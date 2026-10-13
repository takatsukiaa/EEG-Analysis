function results = get_eemd_cca_artifact_removal(eeg_data, Fs, options)
%GET_EEMD_CCA_ARTIFACT_REMOVAL Single-channel EEMD-CCA EMG suppression.
%   RESULTS = GET_EEMD_CCA_ARTIFACT_REMOVAL(EEG_DATA, FS) processes every
%   field of EEG_DATA independently. Ensemble empirical mode decomposition
%   (EEMD) converts each channel into IMF pseudo-channels. CCA between the
%   pseudo-channels and a delayed copy then concentrates temporally
%   decorrelated activity into canonical components.
%
%   An optional synchronized EMG reference can corroborate automatic
%   rejection. This is intentionally supporting evidence: the defining CCA
%   remains between the EEMD components and their delayed copy.
%
%   Important name-value options:
%       EMGReference                 vector or samples-by-EMG matrix ([])
%       WindowDurationSec            default 10
%       OverlapFraction              default 0.5
%       EnsembleSize                 default 10 EMD realizations
%       NoiseStdRatio                default 0.1 of window standard deviation
%       UseComplementaryNoise        default true (+noise/-noise pairs)
%       MaxNumIMF                    default 10
%       DelaySamples                 default 1
%       AutocorrelationThreshold     default 0.9
%       SelectionDirection          "below" for EMG, "above" for the
%                                   high-autocorrelation motion variant
%       RequireEMGConfirmation       default true when EMG is supplied
%       EMGCorrelationThreshold      default 0.10
%       EMGCoherenceThreshold        default 0.10
%       MaxAutoRejectedFraction      default 1 (published threshold rule)
%       ArtifactComponentGain       default 0 (published hard rejection)
%       RejectedComponentsOverride  [] for automatic selection, a numeric
%                                   vector used everywhere, or a cell array
%                                   with one entry per window/channel
%       RandomSeed                   default 0; local deterministic streams
%       StoreDecomposition           store EEMD modes in results (false)
%
%   The implementation follows Sweeney et al., IEEE TBME 2013,
%   doi:10.1109/TBME.2012.2225427. The default low-autocorrelation rule and
%   0.9 threshold follow later single-channel EEG muscle-artifact studies.
%   Complementary noise pairs are based on the public eemd-matlab pattern.
%
%   Fast ictal activity can resemble EMG. Always review the untouched EEG,
%   component diagnostics, and preservation metrics before interpretation.

    arguments
        eeg_data (1,1) struct
        Fs (1,1) double {mustBeFinite,mustBePositive}
        options.EMGReference = []
        options.WindowDurationSec (1,1) double ...
            {mustBeFinite,mustBePositive} = 10
        options.OverlapFraction (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.5
        options.EnsembleSize (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 10
        options.NoiseStdRatio (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.1
        options.UseComplementaryNoise (1,1) logical = true
        options.MaxNumIMF (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 10
        options.SiftRelativeTolerance (1,1) double ...
            {mustBeFinite,mustBePositive} = 0.2
        options.SiftMaxIterations (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 100
        options.DelaySamples (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBePositive} = 1
        options.AutocorrelationThreshold (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.9
        options.SelectionDirection (1,1) string = "below"
        options.RequireEMGConfirmation (1,1) logical = true
        options.EMGCorrelationThreshold (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.10
        options.EMGCoherenceThreshold (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.10
        options.EMGMaxLagSec (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0.050
        options.EMGCoherenceBandHz (1,2) double ...
            {mustBeFinite,mustBeNonnegative} = [30 95]
        options.MaxAutoRejectedFraction (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 1
        options.ArtifactComponentGain (1,1) double ...
            {mustBeFinite,mustBeNonnegative} = 0
        options.RejectedComponentsOverride = []
        options.RandomSeed (1,1) double ...
            {mustBeFinite,mustBeInteger,mustBeNonnegative} = 0
        options.StoreDecomposition (1,1) logical = false
    end

    validate_options(options);
    require_toolboxes();
    [eeg_matrix, channel_names] = eeg_structure_to_matrix(eeg_data);
    [num_samples, num_channels] = size(eeg_matrix);
    emg_reference = resolve_emg_reference(options.EMGReference, num_samples);

    requested_window_samples = max(2, round(options.WindowDurationSec * Fs));
    window_samples = min(num_samples, requested_window_samples);
    minimum_samples = options.DelaySamples + options.MaxNumIMF + 3;
    if window_samples < minimum_samples
        error('get_eemd_cca_artifact_removal:SignalTooShort', ...
            ['Each window must contain at least DelaySamples + MaxNumIMF ' ...
             '+ 3 samples (%d samples for these settings).'], minimum_samples);
    end

    hop_samples = max(1, round(window_samples * (1-options.OverlapFraction)));
    window_starts = make_window_starts(num_samples, window_samples, hop_samples);
    num_windows = numel(window_starts);
    overrides = resolve_overrides(options.RejectedComponentsOverride, ...
        num_windows, num_channels);

    cleaned_accumulator = zeros(size(eeg_matrix));
    weight_accumulator = zeros(num_samples, 1);
    window_ends = zeros(num_windows, 1);
    component_autocorrelations = cell(num_windows, num_channels);
    component_features = cell(num_windows, num_channels);
    component_gains = cell(num_windows, num_channels);
    rejected_components = cell(num_windows, num_channels);
    component_unmixing = cell(num_windows, num_channels);
    component_mixing = cell(num_windows, num_channels);
    eemd_trial_imf_counts = cell(num_windows, num_channels);
    eemd_component_count = zeros(num_windows, num_channels);
    effective_component_rank = zeros(num_windows, num_channels);
    eemd_reconstruction_error = zeros(num_windows, num_channels);
    cca_reconstruction_error = zeros(num_windows, num_channels);
    skipped_reason = strings(num_windows, num_channels);
    if options.StoreDecomposition
        eemd_components = cell(num_windows, num_channels);
    else
        eemd_components = cell(0, 0);
    end

    for window_index = 1:num_windows
        first_sample = window_starts(window_index);
        last_sample = min(first_sample+window_samples-1, num_samples);
        indices = first_sample:last_sample;
        window_ends(window_index) = last_sample;
        if isempty(emg_reference)
            emg_window = [];
        else
            emg_window = emg_reference(indices, :);
        end
        window_weight = overlap_weight(numel(indices));

        for channel_index = 1:num_channels
            signal = eeg_matrix(indices, channel_index);
            stream_seed = mod(options.RandomSeed + ...
                (window_index-1)*num_channels + channel_index-1, 2^32-1);
            [cleaned, diagnostics] = clean_eemd_cca_window( ...
                signal, emg_window, Fs, stream_seed, ...
                overrides{window_index, channel_index}, options);

            cleaned_accumulator(indices, channel_index) = ...
                cleaned_accumulator(indices, channel_index) + ...
                cleaned .* window_weight;
            component_autocorrelations{window_index, channel_index} = ...
                diagnostics.component_autocorrelations;
            component_features{window_index, channel_index} = ...
                diagnostics.component_features;
            component_gains{window_index, channel_index} = ...
                diagnostics.component_gains;
            rejected_components{window_index, channel_index} = ...
                diagnostics.rejected_components;
            component_unmixing{window_index, channel_index} = ...
                diagnostics.component_unmixing;
            component_mixing{window_index, channel_index} = ...
                diagnostics.component_mixing;
            eemd_trial_imf_counts{window_index, channel_index} = ...
                diagnostics.eemd_trial_imf_counts;
            eemd_component_count(window_index, channel_index) = ...
                diagnostics.eemd_component_count;
            effective_component_rank(window_index, channel_index) = ...
                diagnostics.effective_component_rank;
            eemd_reconstruction_error(window_index, channel_index) = ...
                diagnostics.eemd_reconstruction_error;
            cca_reconstruction_error(window_index, channel_index) = ...
                diagnostics.cca_reconstruction_error;
            skipped_reason(window_index, channel_index) = ...
                diagnostics.skipped_reason;
            if options.StoreDecomposition
                eemd_components{window_index, channel_index} = ...
                    diagnostics.eemd_components;
            end
        end
        weight_accumulator(indices) = weight_accumulator(indices) + window_weight;
    end

    if any(weight_accumulator <= 0)
        error('get_eemd_cca_artifact_removal:InvalidWindowCoverage', ...
            'The overlap-add windows did not cover every output sample.');
    end
    cleaned_matrix = cleaned_accumulator ./ weight_accumulator;
    artifact_matrix = eeg_matrix-cleaned_matrix;
    maximum_emg_lag = round(options.EMGMaxLagSec*Fs);

    results = struct();
    results.cleaned_eeg = matrix_to_eeg_structure(cleaned_matrix, channel_names);
    results.artifact_estimate = matrix_to_eeg_structure(artifact_matrix, channel_names);
    results.channel_names = channel_names;
    results.window_start_samples = window_starts(:);
    results.window_end_samples = window_ends;
    results.window_start_sec = (window_starts(:)-1)/Fs;
    results.window_end_sec = window_ends/Fs;
    results.component_autocorrelations = component_autocorrelations;
    results.component_features = component_features;
    results.component_gains = component_gains;
    results.rejected_components = rejected_components;
    results.component_unmixing = component_unmixing;
    results.component_mixing = component_mixing;
    results.eemd_trial_imf_counts = eemd_trial_imf_counts;
    results.eemd_component_count = eemd_component_count;
    results.effective_component_rank = effective_component_rank;
    results.eemd_reconstruction_error = eemd_reconstruction_error;
    results.cca_reconstruction_error = cca_reconstruction_error;
    results.skipped_reason = skipped_reason;
    if options.StoreDecomposition
        results.eemd_components = eemd_components;
    end
    results.metrics = calculate_metrics(eeg_matrix, cleaned_matrix, ...
        artifact_matrix, emg_reference, maximum_emg_lag, channel_names);
    results.settings = struct( ...
        'method', 'Windowed single-channel EEMD-CCA with optional EMG confirmation', ...
        'window_duration_sec', options.WindowDurationSec, ...
        'window_samples', window_samples, ...
        'overlap_fraction', options.OverlapFraction, ...
        'hop_samples', hop_samples, ...
        'ensemble_size', options.EnsembleSize, ...
        'noise_standard_deviation_ratio', options.NoiseStdRatio, ...
        'use_complementary_noise', options.UseComplementaryNoise, ...
        'maximum_imf_count', options.MaxNumIMF, ...
        'sift_relative_tolerance', options.SiftRelativeTolerance, ...
        'sift_maximum_iterations', options.SiftMaxIterations, ...
        'delay_samples', options.DelaySamples, ...
        'delay_sec', options.DelaySamples/Fs, ...
        'autocorrelation_threshold', options.AutocorrelationThreshold, ...
        'selection_direction', char(options.SelectionDirection), ...
        'has_emg_reference', ~isempty(emg_reference), ...
        'require_emg_confirmation', options.RequireEMGConfirmation, ...
        'emg_correlation_threshold', options.EMGCorrelationThreshold, ...
        'emg_coherence_threshold', options.EMGCoherenceThreshold, ...
        'emg_maximum_lag_sec', options.EMGMaxLagSec, ...
        'emg_coherence_band_hz', options.EMGCoherenceBandHz, ...
        'max_auto_rejected_fraction', options.MaxAutoRejectedFraction, ...
        'artifact_component_gain', options.ArtifactComponentGain, ...
        'automatic_selection', isempty(options.RejectedComponentsOverride), ...
        'random_seed', options.RandomSeed, ...
        'store_decomposition', options.StoreDecomposition);
    results.notes = [ ...
        "EEMD-CCA is applied independently to each EEG channel.", ...
        "EMG corroborates selection but is not an input to the CCA decomposition.", ...
        "Fast ictal EEG can resemble EMG; compare with the untouched recording."];
end

function validate_options(options)
    if options.OverlapFraction >= 1
        error('get_eemd_cca_artifact_removal:InvalidOverlap', ...
            'OverlapFraction must be in [0,1).');
    end
    bounded = {'AutocorrelationThreshold','EMGCorrelationThreshold', ...
        'EMGCoherenceThreshold','MaxAutoRejectedFraction', ...
        'ArtifactComponentGain'};
    for index = 1:numel(bounded)
        if options.(bounded{index}) > 1
            error('get_eemd_cca_artifact_removal:InvalidOption', ...
                '%s must be in [0,1].', bounded{index});
        end
    end
    direction = lower(options.SelectionDirection);
    if direction ~= "below" && direction ~= "above"
        error('get_eemd_cca_artifact_removal:InvalidSelectionDirection', ...
            'SelectionDirection must be "below" or "above".');
    end
    band = options.EMGCoherenceBandHz;
    if band(1) >= band(2)
        error('get_eemd_cca_artifact_removal:InvalidEMGBand', ...
            'EMGCoherenceBandHz must be [low high] with low < high.');
    end
end

function require_toolboxes()
    if exist('emd', 'file') ~= 2
        error('get_eemd_cca_artifact_removal:MissingSignalToolbox', ...
            'Signal Processing Toolbox is required; emd was not found.');
    end
    if exist('canoncorr', 'file') ~= 2
        error('get_eemd_cca_artifact_removal:MissingStatisticsToolbox', ...
            ['Statistics and Machine Learning Toolbox is required; ' ...
             'canoncorr was not found.']);
    end
end

function [cleaned, diagnostics] = clean_eemd_cca_window( ...
        signal, emg, Fs, stream_seed, override, options)
    cleaned = signal;
    diagnostics = empty_window_diagnostics();
    [components, trial_counts] = ensemble_emd(signal, stream_seed, options);
    diagnostics.eemd_components = components;
    diagnostics.eemd_trial_imf_counts = trial_counts;
    diagnostics.eemd_component_count = size(components, 2);
    diagnostics.eemd_reconstruction_error = norm(sum(components,2)-signal) / ...
        max(norm(signal), eps);

    if size(components, 2) < 2
        diagnostics.skipped_reason = "EEMD returned fewer than two components";
        return;
    end

    component_center = mean(components, 1);
    centered = components-component_center;
    current = centered(options.DelaySamples+1:end, :);
    [basis, component_rank] = right_singular_basis(current);
    diagnostics.effective_component_rank = component_rank;
    if component_rank < 1
        diagnostics.skipped_reason = "EEMD components have zero numerical rank";
        return;
    end
    reduced = centered*basis;
    try
        [cca_coefficients, ~, correlations] = canoncorr( ...
            reduced(options.DelaySamples+1:end,:), ...
            reduced(1:end-options.DelaySamples,:));
    catch exception
        diagnostics.skipped_reason = "CCA failed: " + string(exception.message);
        return;
    end
    correlations = double(correlations(:));
    canonical_count = numel(correlations);
    if canonical_count < 1
        diagnostics.skipped_reason = "CCA returned no canonical components";
        return;
    end
    if size(cca_coefficients,2) < component_rank
        transform = [cca_coefficients null(cca_coefficients.')];
    else
        transform = cca_coefficients;
    end
    if size(transform,2) ~= component_rank || rank(transform) ~= component_rank
        diagnostics.skipped_reason = "CCA transform is numerically singular";
        return;
    end

    autocorrelations = nan(component_rank,1);
    autocorrelations(1:canonical_count) = correlations;
    unmixing = basis*transform;
    mixing = pinv(unmixing);
    scores = centered*unmixing;
    rank_residual = centered-scores*mixing;
    no_rejection_components = scores*mixing+rank_residual+component_center;
    no_rejection = sum(no_rejection_components,2);
    diagnostics.cca_reconstruction_error = norm(no_rejection-signal) / ...
        max(norm(signal),eps);

    [features, emg_supported] = component_reference_features( ...
        scores, autocorrelations, emg, Fs, options);
    if isempty(override)
        candidates = automatic_candidates(autocorrelations, ...
            emg_supported, ~isempty(emg), options);
    else
        candidates = validate_override(override, component_rank);
    end
    gains = ones(component_rank,1);
    gains(candidates) = options.ArtifactComponentGain;
    features.Selected = ismember(features.Component, candidates(:));
    features.Gain = gains;

    cleaned_components = (scores.*gains.')*mixing + ...
        rank_residual+component_center;
    cleaned = sum(cleaned_components,2);
    if any(~isfinite(cleaned))
        error('get_eemd_cca_artifact_removal:NonfiniteReconstruction', ...
            'EEMD-CCA reconstruction produced NaN or Inf.');
    end
    diagnostics.component_autocorrelations = autocorrelations;
    diagnostics.component_features = features;
    diagnostics.component_gains = gains;
    diagnostics.rejected_components = candidates(:).';
    diagnostics.component_unmixing = unmixing;
    diagnostics.component_mixing = mixing;
    diagnostics.skipped_reason = "";
end

function [components, trial_counts] = ensemble_emd(signal, stream_seed, options)
    signal = double(signal(:));
    signal_mean = mean(signal);
    signal_scale = std(signal,0);
    if signal_scale <= eps(max(1,abs(signal_mean)))
        components = signal;
        trial_counts = zeros(options.EnsembleSize,1);
        return;
    end
    normalized = (signal-signal_mean)/signal_scale;
    num_samples = numel(signal);
    imf_accumulator = zeros(num_samples, options.MaxNumIMF);
    trial_counts = zeros(options.EnsembleSize,1);
    stream = RandStream('mt19937ar','Seed',double(stream_seed));
    paired_noise = [];

    for trial = 1:options.EnsembleSize
        if options.UseComplementaryNoise && mod(trial,2)==0
            noise = -paired_noise;
        else
            noise = options.NoiseStdRatio*randn(stream,num_samples,1);
            paired_noise = noise;
        end
        noisy = normalized+noise;
        [imfs, ~] = emd(noisy, ...
            'MaxNumIMF',options.MaxNumIMF, ...
            'SiftRelativeTolerance',options.SiftRelativeTolerance, ...
            'SiftMaxIterations',options.SiftMaxIterations, ...
            'Display',0);
        if isempty(imfs)
            continue;
        end
        if size(imfs,1) ~= num_samples && size(imfs,2) == num_samples
            imfs = imfs.';
        end
        count = min(size(imfs,2), options.MaxNumIMF);
        trial_counts(trial) = count;
        imf_accumulator(:,1:count) = imf_accumulator(:,1:count) + ...
            imfs(:,1:count);
    end

    imf_count = max(trial_counts);
    if imf_count == 0
        components = signal;
        return;
    end
    averaged_imfs = signal_scale * ...
        imf_accumulator(:,1:imf_count)/options.EnsembleSize;
    residual = signal-sum(averaged_imfs,2);
    residual = residual+0*signal_mean; % documents that the DC term is retained here
    components = [averaged_imfs residual];
end

function [features, supported] = component_reference_features( ...
        scores, autocorrelations, emg, Fs, options)
    count = size(scores,2);
    maximum_lag = round(options.EMGMaxLagSec*Fs);
    emg_correlation = nan(count,1);
    emg_coherence = nan(count,1);
    if ~isempty(emg)
        for component = 1:count
            emg_correlation(component) = max_reference_correlation( ...
                scores(:,component), emg, maximum_lag);
            emg_coherence(component) = max_reference_coherence( ...
                scores(:,component), emg, Fs, options.EMGCoherenceBandHz);
        end
    end
    supported = emg_correlation >= options.EMGCorrelationThreshold | ...
        emg_coherence >= options.EMGCoherenceThreshold;
    supported(isnan(supported)) = false;
    features = table((1:count)', autocorrelations, emg_correlation, ...
        emg_coherence, supported, ...
        'VariableNames', {'Component','Autocorrelation', ...
        'MaxAbsEMGCorrelation','MuscleBandEMGCoherence','EMGSupported'});
end

function candidates = automatic_candidates(autocorrelations, ...
        emg_supported, has_emg, options)
    valid = isfinite(autocorrelations);
    if lower(options.SelectionDirection) == "below"
        mask = valid & autocorrelations < options.AutocorrelationThreshold;
        [~, order] = sort(autocorrelations, 'ascend', 'MissingPlacement','last');
    else
        mask = valid & autocorrelations > options.AutocorrelationThreshold;
        [~, order] = sort(autocorrelations, 'descend', 'MissingPlacement','last');
    end
    if has_emg && options.RequireEMGConfirmation
        mask = mask & emg_supported;
    end
    ordered = order(mask(order));
    if options.MaxAutoRejectedFraction == 0 || isempty(ordered)
        candidates = zeros(1,0);
        return;
    end
    maximum = max(1, floor(options.MaxAutoRejectedFraction*numel(autocorrelations)));
    candidates = sort(ordered(1:min(maximum,numel(ordered)))).';
end

function value = max_reference_correlation(signal, references, maximum_lag)
    value = 0;
    for reference_index = 1:size(references,2)
        reference = references(:,reference_index);
        for lag = -maximum_lag:maximum_lag
            if lag >= 0
                first = signal(1+lag:end);
                second = reference(1:end-lag);
            else
                first = signal(1:end+lag);
                second = reference(1-lag:end);
            end
            value = max(value, absolute_correlation(first,second));
        end
    end
end

function value = max_reference_coherence(signal, references, Fs, band)
    upper = min(band(2), Fs/2*(1-1e-6));
    lower = band(1);
    if upper <= lower || numel(signal) < 32
        value = NaN;
        return;
    end
    window_length = min(numel(signal), max(32,round(2*Fs)));
    overlap = floor(window_length/2);
    nfft = max(256,2^nextpow2(window_length));
    value = 0;
    for reference_index = 1:size(references,2)
        [coherence, frequencies] = mscohere(signal, ...
            references(:,reference_index), hamming(window_length), ...
            overlap, nfft, Fs);
        mask = frequencies >= lower & frequencies <= upper;
        if any(mask)
            value = max(value, mean(coherence(mask)));
        end
    end
end

function reference = resolve_emg_reference(requested, num_samples)
    if isempty(requested)
        reference = [];
        return;
    end
    if isstruct(requested)
        names = fieldnames(requested);
        reference = zeros(num_samples,numel(names));
        for index = 1:numel(names)
            signal = requested.(names{index});
            if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                    numel(signal) ~= num_samples || any(~isfinite(signal(:)))
                error('get_eemd_cca_artifact_removal:InvalidEMGReference', ...
                    'Every EMG structure field must be a finite %d-sample vector.', ...
                    num_samples);
            end
            reference(:,index) = double(signal(:));
        end
    elseif isnumeric(requested) && isreal(requested) && ...
            ismatrix(requested) && all(isfinite(requested(:)))
        if isvector(requested) && numel(requested) == num_samples
            reference = double(requested(:));
        elseif size(requested,1) == num_samples
            reference = double(requested);
        elseif size(requested,2) == num_samples
            reference = double(requested.');
        else
            error('get_eemd_cca_artifact_removal:InvalidEMGReference', ...
                'EMGReference must have %d samples.',num_samples);
        end
    else
        error('get_eemd_cca_artifact_removal:InvalidEMGReference', ...
            'EMGReference must be a finite real vector, matrix, or structure.');
    end
    reference = reference-mean(reference,1);
end

function overrides = resolve_overrides(requested, num_windows, num_channels)
    overrides = repmat({[]},num_windows,num_channels);
    if isempty(requested)
        return;
    end
    if isnumeric(requested)
        overrides(:) = {validate_override_shape(requested)};
    elseif iscell(requested) && numel(requested) == num_windows*num_channels
        requested = reshape(requested,num_windows,num_channels);
        for index = 1:numel(requested)
            overrides{index} = validate_override_shape(requested{index});
        end
    else
        error('get_eemd_cca_artifact_removal:InvalidOverride', ...
            ['RejectedComponentsOverride must be empty, a numeric vector, ' ...
             'or a cell array with one entry per window/channel.']);
    end
end

function indices = validate_override_shape(indices)
    if isempty(indices)
        indices = zeros(1,0);
    elseif ~isnumeric(indices) || ~isreal(indices) || ~isvector(indices) || ...
            any(~isfinite(indices)) || any(indices ~= round(indices)) || ...
            any(indices < 1)
        error('get_eemd_cca_artifact_removal:InvalidOverride', ...
            'Rejected component indices must be positive finite integers.');
    else
        indices = unique(double(indices(:).'),'stable');
    end
end

function indices = validate_override(indices, maximum)
    indices = validate_override_shape(indices);
    if any(indices > maximum)
        error('get_eemd_cca_artifact_removal:InvalidOverrideForWindow', ...
            'A rejected component index exceeds this window rank (%d).',maximum);
    end
end

function [basis, effective_rank] = right_singular_basis(matrix)
    [~, singular_values, right_vectors] = svd(matrix,'econ');
    values = diag(singular_values);
    if isempty(values) || max(values) <= eps
        basis = zeros(size(matrix,2),0);
        effective_rank = 0;
        return;
    end
    tolerance = max(size(matrix))*eps(max(values));
    effective_rank = sum(values > tolerance);
    basis = right_vectors(:,1:effective_rank);
end

function diagnostics = empty_window_diagnostics()
    diagnostics = struct( ...
        'component_autocorrelations',zeros(0,1), ...
        'component_features',table(), ...
        'component_gains',zeros(0,1), ...
        'rejected_components',zeros(1,0), ...
        'component_unmixing',zeros(0,0), ...
        'component_mixing',zeros(0,0), ...
        'eemd_trial_imf_counts',zeros(0,1), ...
        'eemd_component_count',0, ...
        'effective_component_rank',0, ...
        'eemd_reconstruction_error',0, ...
        'cca_reconstruction_error',0, ...
        'eemd_components',zeros(0,0), ...
        'skipped_reason',"");
end

function metrics = calculate_metrics(original, cleaned, artifact, ...
        emg, maximum_lag, channel_names)
    count = numel(channel_names);
    original_rms = zeros(count,1);
    cleaned_rms = zeros(count,1);
    artifact_rms = zeros(count,1);
    artifact_ratio = zeros(count,1);
    preservation_correlation = zeros(count,1);
    original_emg_correlation = nan(count,1);
    artifact_emg_correlation = nan(count,1);
    for channel = 1:count
        before = original(:,channel)-mean(original(:,channel));
        after = cleaned(:,channel)-mean(cleaned(:,channel));
        removed = artifact(:,channel);
        original_rms(channel) = root_mean_square(before);
        cleaned_rms(channel) = root_mean_square(after);
        artifact_rms(channel) = root_mean_square(removed);
        artifact_ratio(channel) = safe_ratio(artifact_rms(channel),original_rms(channel));
        preservation_correlation(channel) = absolute_correlation(before,after);
        if ~isempty(emg)
            original_emg_correlation(channel) = ...
                max_reference_correlation(before,emg,maximum_lag);
            artifact_emg_correlation(channel) = ...
                max_reference_correlation(removed,emg,maximum_lag);
        end
    end
    rms_change_percent = 100*(cleaned_rms-original_rms)./max(original_rms,eps);
    metrics = table(channel_names,original_rms,cleaned_rms,artifact_rms, ...
        artifact_ratio,preservation_correlation,original_emg_correlation, ...
        artifact_emg_correlation,rms_change_percent, ...
        'VariableNames',{'Channel','OriginalRMS','CleanedRMS','ArtifactRMS', ...
        'ArtifactToSignalRMSRatio','OriginalToCleanedCorrelation', ...
        'OriginalMaxAbsEMGCorrelation','ArtifactMaxAbsEMGCorrelation', ...
        'RMSChangePercent'});
end

function [matrix, channel_names] = eeg_structure_to_matrix(eeg_data)
    channel_names = fieldnames(eeg_data);
    if isempty(channel_names)
        error('get_eemd_cca_artifact_removal:NoEEGChannels', ...
            'eeg_data must contain at least one channel.');
    end
    num_samples = numel(eeg_data.(channel_names{1}));
    matrix = zeros(num_samples,numel(channel_names));
    for channel = 1:numel(channel_names)
        signal = eeg_data.(channel_names{channel});
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('get_eemd_cca_artifact_removal:InvalidEEGChannel', ...
                'EEG channel "%s" must be a finite real vector.', ...
                channel_names{channel});
        end
        if numel(signal) ~= num_samples
            error('get_eemd_cca_artifact_removal:UnequalChannelLengths', ...
                'EEG channel "%s" has an inconsistent length.', ...
                channel_names{channel});
        end
        matrix(:,channel) = double(signal(:));
    end
end

function structure = matrix_to_eeg_structure(matrix, channel_names)
    structure = struct();
    for channel = 1:numel(channel_names)
        structure.(channel_names{channel}) = matrix(:,channel);
    end
end

function starts = make_window_starts(num_samples,window_samples,hop_samples)
    if num_samples <= window_samples
        starts = 1;
        return;
    end
    last_start = num_samples-window_samples+1;
    starts = 1:hop_samples:last_start;
    if starts(end) ~= last_start
        starts(end+1) = last_start;
    end
end

function weight = overlap_weight(num_samples)
    phase = pi*((1:num_samples)'-0.5)/num_samples;
    weight = sin(phase).^2;
end

function value = absolute_correlation(first,second)
    first = double(first(:))-mean(first);
    second = double(second(:))-mean(second);
    denominator = norm(first)*norm(second);
    if denominator <= eps
        value = 0;
    else
        value = abs((first'*second)/denominator);
    end
end

function value = root_mean_square(signal)
    value = sqrt(mean(double(signal(:)).^2));
end

function value = safe_ratio(numerator,denominator)
    if denominator <= eps(max(1,abs(numerator)))
        if abs(numerator) <= eps
            value = 0;
        else
            value = Inf;
        end
    else
        value = numerator/denominator;
    end
end
