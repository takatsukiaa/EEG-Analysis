function results = get_eeg_bss_cca_artifact_removal(eeg_data, Fs, options)
%GET_EEG_BSS_CCA_ARTIFACT_REMOVAL Remove low-autocorrelation EEG sources.
%   RESULTS = GET_EEG_BSS_CCA_ARTIFACT_REMOVAL(EEG_DATA, FS) applies
%   windowed blind-source-separation CCA between multichannel EEG and a
%   delayed copy of the same EEG. Low-autocorrelation canonical components
%   are treated as muscle-like artifacts and excluded from reconstruction.
%
%   Name-value options:
%       WindowDurationSec          - CCA window duration (default 10)
%       OverlapFraction            - window overlap in [0,1) (default 0.5)
%       DelaySamples               - delayed-EEG offset (default 1)
%       AutocorrelationThreshold   - reject components below this value
%                                    (default 0.80)
%       MaxAutoRejectedFraction    - automatic per-window cap (default 0.25)
%       RejectedComponentsOverride - [] for automatic selection, a numeric
%                                    vector applied to every window, or one
%                                    vector per window in a cell array
%
%   Low-autocorrelation ictal activity can be mistaken for muscle. Compare
%   every result with the untouched recording before interpretation.

    arguments
        eeg_data (1, 1) struct
        Fs (1, 1) double {mustBeFinite, mustBePositive}
        options.WindowDurationSec (1, 1) double ...
            {mustBeFinite, mustBePositive} = 10
        options.OverlapFraction (1, 1) double ...
            {mustBeFinite, mustBeNonnegative} = 0.5
        options.DelaySamples (1, 1) double ...
            {mustBeFinite, mustBePositive, mustBeInteger} = 1
        options.AutocorrelationThreshold (1, 1) double ...
            {mustBeFinite, mustBeNonnegative} = 0.80
        options.MaxAutoRejectedFraction (1, 1) double ...
            {mustBeFinite, mustBeNonnegative} = 0.25
        options.RejectedComponentsOverride = []
    end

    if options.OverlapFraction >= 1
        error('get_eeg_bss_cca_artifact_removal:InvalidOverlap', ...
            'OverlapFraction must be greater than or equal to 0 and less than 1.');
    end
    if options.AutocorrelationThreshold > 1
        error('get_eeg_bss_cca_artifact_removal:InvalidThreshold', ...
            'AutocorrelationThreshold must be in the interval [0, 1].');
    end
    if options.MaxAutoRejectedFraction > 1
        error('get_eeg_bss_cca_artifact_removal:InvalidRejectedFraction', ...
            'MaxAutoRejectedFraction must be in the interval [0, 1].');
    end
    if exist('canoncorr', 'file') ~= 2
        error('get_eeg_bss_cca_artifact_removal:MissingStatisticsToolbox', ...
            ['Statistics and Machine Learning Toolbox is required. ' ...
             'The function canoncorr was not found on the MATLAB path.']);
    end

    [eeg_matrix, channel_names] = eeg_structure_to_matrix(eeg_data);
    num_samples = size(eeg_matrix, 1);
    num_channels = size(eeg_matrix, 2);
    if num_channels < 2
        error('get_eeg_bss_cca_artifact_removal:NotEnoughEEGChannels', ...
            'At least two EEG channels are required for multichannel BSS-CCA.');
    end

    requested_window_samples = max(2, round(options.WindowDurationSec * Fs));
    window_samples = min(num_samples, requested_window_samples);
    if window_samples <= options.DelaySamples + num_channels
        error('get_eeg_bss_cca_artifact_removal:SignalTooShort', ...
            ['Each CCA window must contain more than the delay (%d samples) ' ...
             'plus the number of EEG channels (%d).'], ...
            options.DelaySamples, num_channels);
    end

    hop_samples = max(1, round(window_samples * ...
        (1 - options.OverlapFraction)));
    window_starts = make_window_starts(num_samples, window_samples, hop_samples);
    num_windows = numel(window_starts);
    [override_by_window, automatic_selection] = resolve_overrides( ...
        options.RejectedComponentsOverride, num_windows, num_channels);

    cleaned_accumulator = zeros(size(eeg_matrix));
    weight_accumulator = zeros(num_samples, 1);
    window_ends = zeros(num_windows, 1);
    canonical_autocorrelations = cell(num_windows, 1);
    rejected_components = cell(num_windows, 1);
    component_unmixing = cell(num_windows, 1);
    component_mixing = cell(num_windows, 1);
    effective_eeg_rank = zeros(num_windows, 1);
    reconstruction_error = zeros(num_windows, 1);
    skipped_reason = strings(num_windows, 1);

    for window_index = 1:num_windows
        first_sample = window_starts(window_index);
        last_sample = min(first_sample + window_samples - 1, num_samples);
        window_ends(window_index) = last_sample;
        indices = first_sample:last_sample;
        eeg_window = eeg_matrix(indices, :);

        [cleaned_window, diagnostics] = clean_bss_cca_window( ...
            eeg_window, options.DelaySamples, ...
            options.AutocorrelationThreshold, ...
            options.MaxAutoRejectedFraction, automatic_selection, ...
            override_by_window{window_index});

        window_weight = overlap_weight(numel(indices));
        cleaned_accumulator(indices, :) = cleaned_accumulator(indices, :) + ...
            cleaned_window .* window_weight;
        weight_accumulator(indices) = weight_accumulator(indices) + window_weight;

        canonical_autocorrelations{window_index} = ...
            diagnostics.canonical_autocorrelations;
        rejected_components{window_index} = diagnostics.rejected_components;
        component_unmixing{window_index} = diagnostics.component_unmixing;
        component_mixing{window_index} = diagnostics.component_mixing;
        effective_eeg_rank(window_index) = diagnostics.effective_eeg_rank;
        reconstruction_error(window_index) = diagnostics.reconstruction_error;
        skipped_reason(window_index) = diagnostics.skipped_reason;
    end

    if any(weight_accumulator <= 0)
        error('get_eeg_bss_cca_artifact_removal:InvalidWindowCoverage', ...
            'The overlap-add windows did not cover every output sample.');
    end
    cleaned_matrix = cleaned_accumulator ./ weight_accumulator;
    artifact_matrix = eeg_matrix - cleaned_matrix;

    results = struct();
    results.cleaned_eeg = matrix_to_eeg_structure(cleaned_matrix, channel_names);
    results.artifact_estimate = ...
        matrix_to_eeg_structure(artifact_matrix, channel_names);
    results.channel_names = channel_names;
    results.window_start_samples = window_starts(:);
    results.window_end_samples = window_ends;
    results.window_start_sec = (window_starts(:) - 1) / Fs;
    results.window_end_sec = window_ends / Fs;
    results.canonical_autocorrelations = canonical_autocorrelations;
    results.rejected_components = rejected_components;
    results.component_unmixing = component_unmixing;
    results.component_mixing = component_mixing;
    results.effective_eeg_rank = effective_eeg_rank;
    results.reconstruction_error = reconstruction_error;
    results.skipped_reason = skipped_reason;
    results.metrics = calculate_metrics( ...
        eeg_matrix, cleaned_matrix, artifact_matrix, channel_names, ...
        options.DelaySamples);
    results.settings = struct( ...
        'method', 'EEG delayed-copy BSS-CCA', ...
        'window_duration_sec', options.WindowDurationSec, ...
        'window_samples', window_samples, ...
        'overlap_fraction', options.OverlapFraction, ...
        'hop_samples', hop_samples, ...
        'delay_samples', options.DelaySamples, ...
        'delay_sec', options.DelaySamples / Fs, ...
        'autocorrelation_threshold', options.AutocorrelationThreshold, ...
        'max_auto_rejected_fraction', options.MaxAutoRejectedFraction, ...
        'automatic_selection', automatic_selection, ...
        'rejected_components_override', ...
            options.RejectedComponentsOverride);
    results.notes = [ ...
        "Low-autocorrelation ictal activity can be mistaken for muscle. " ...
        "Compare cleaned signals with the untouched recording before interpretation."];
end

function [matrix, channel_names] = eeg_structure_to_matrix(eeg_data)
    channel_names = fieldnames(eeg_data);
    if isempty(channel_names)
        error('get_eeg_bss_cca_artifact_removal:NoEEGChannels', ...
            'eeg_data must contain at least one channel.');
    end
    first_signal = eeg_data.(channel_names{1});
    if ~isnumeric(first_signal) || ~isreal(first_signal) || ...
            ~isvector(first_signal) || isempty(first_signal)
        error('get_eeg_bss_cca_artifact_removal:InvalidEEGChannel', ...
            'Every EEG channel must be a nonempty finite real numeric vector.');
    end
    num_samples = numel(first_signal);
    matrix = zeros(num_samples, numel(channel_names));
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = eeg_data.(channel_name);
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('get_eeg_bss_cca_artifact_removal:InvalidEEGChannel', ...
                'EEG channel "%s" must be a finite real numeric vector.', ...
                channel_name);
        end
        if numel(signal) ~= num_samples
            error('get_eeg_bss_cca_artifact_removal:UnequalChannelLengths', ...
                'EEG channel "%s" has an inconsistent length.', channel_name);
        end
        matrix(:, channel_index) = double(signal(:));
    end
end

function structure = matrix_to_eeg_structure(matrix, channel_names)
    structure = struct();
    for channel_index = 1:numel(channel_names)
        structure.(channel_names{channel_index}) = matrix(:, channel_index);
    end
end

function starts = make_window_starts(num_samples, window_samples, hop_samples)
    if num_samples <= window_samples
        starts = 1;
        return;
    end
    last_start = num_samples - window_samples + 1;
    starts = 1:hop_samples:last_start;
    if starts(end) ~= last_start
        starts(end + 1) = last_start;
    end
end

function [overrides, automatic] = resolve_overrides( ...
        requested, num_windows, num_channels)
    automatic = isempty(requested);
    overrides = repmat({[]}, num_windows, 1);
    if automatic
        return;
    end
    if isnumeric(requested)
        validated = validate_component_indices(requested, num_channels);
        overrides(:) = {validated};
    elseif iscell(requested)
        if numel(requested) ~= num_windows
            error('get_eeg_bss_cca_artifact_removal:InvalidOverride', ...
                ['A cell RejectedComponentsOverride must contain one entry ' ...
                 'for each of the %d windows.'], num_windows);
        end
        for window_index = 1:num_windows
            overrides{window_index} = validate_component_indices( ...
                requested{window_index}, num_channels);
        end
    else
        error('get_eeg_bss_cca_artifact_removal:InvalidOverride', ...
            ['RejectedComponentsOverride must be empty, a numeric vector, ' ...
             'or a cell array of numeric vectors.']);
    end
end

function indices = validate_component_indices(indices, maximum_index)
    if isempty(indices)
        indices = zeros(1, 0);
        return;
    end
    if ~isnumeric(indices) || ~isreal(indices) || ~isvector(indices) || ...
            any(~isfinite(indices)) || any(indices ~= round(indices)) || ...
            any(indices < 1) || any(indices > maximum_index)
        error('get_eeg_bss_cca_artifact_removal:InvalidOverride', ...
            'Rejected component indices must be unique integers from 1 to %d.', ...
            maximum_index);
    end
    indices = unique(double(indices(:).'), 'stable');
end

function [cleaned, diagnostics] = clean_bss_cca_window( ...
        eeg, delay_samples, threshold, maximum_fraction, ...
        automatic_selection, override)
    [num_samples, num_channels] = size(eeg);
    diagnostics = empty_window_diagnostics();
    cleaned = eeg;

    current_eeg = eeg(delay_samples + 1:end, :);
    eeg_center = mean(current_eeg, 1);
    centered_eeg = eeg - eeg_center;
    [eeg_basis, eeg_rank] = right_singular_basis( ...
        centered_eeg(delay_samples + 1:end, :));
    diagnostics.effective_eeg_rank = eeg_rank;
    if eeg_rank == 0
        diagnostics.skipped_reason = "EEG window has zero numerical rank";
        return;
    end

    reduced_eeg = centered_eeg * eeg_basis;
    [canonical_coefficients, ~, correlations] = canoncorr( ...
        reduced_eeg(delay_samples + 1:end, :), ...
        reduced_eeg(1:end - delay_samples, :));
    correlations = double(correlations(:));
    canonical_count = numel(correlations);
    if canonical_count == 0
        diagnostics.skipped_reason = "CCA returned no canonical components";
        return;
    end

    if size(canonical_coefficients, 2) < eeg_rank
        complement = null(canonical_coefficients.');
        full_transform = [canonical_coefficients complement];
    else
        full_transform = canonical_coefficients;
    end
    if size(full_transform, 2) ~= eeg_rank || ...
            rank(full_transform) ~= eeg_rank
        diagnostics.skipped_reason = "CCA transform is numerically singular";
        return;
    end

    full_correlations = zeros(eeg_rank, 1);
    full_correlations(1:canonical_count) = correlations;
    if automatic_selection
        candidates = find(correlations < threshold).';
        if maximum_fraction == 0 || isempty(candidates)
            rejected = zeros(1, 0);
        else
            maximum_rejected = max(1, floor(maximum_fraction * eeg_rank));
            [~, candidate_order] = sort(correlations(candidates), 'ascend');
            candidates = candidates(candidate_order);
            rejected = candidates(1:min(numel(candidates), maximum_rejected));
            rejected = sort(rejected);
        end
    else
        if any(override > eeg_rank)
            error('get_eeg_bss_cca_artifact_removal:InvalidOverrideForWindow', ...
                ['A rejected component index exceeds the effective EEG ' ...
                 'rank (%d) of this window.'], eeg_rank);
        end
        rejected = override;
    end

    scores = reduced_eeg * full_transform;
    inverse_transform = pinv(full_transform);
    no_rejection_reduced = scores * inverse_transform;
    rank_residual = centered_eeg - reduced_eeg * eeg_basis.';
    reconstructed = no_rejection_reduced * eeg_basis.' + ...
        rank_residual + eeg_center;
    reconstruction_denominator = max(norm(eeg, 'fro'), eps);
    diagnostics.reconstruction_error = ...
        norm(reconstructed - eeg, 'fro') / reconstruction_denominator;

    scores(:, rejected) = 0;
    cleaned_reduced = scores * inverse_transform;
    cleaned = cleaned_reduced * eeg_basis.' + rank_residual + eeg_center;
    if any(~isfinite(cleaned), 'all')
        error('get_eeg_bss_cca_artifact_removal:NonfiniteReconstruction', ...
            'BSS-CCA reconstruction produced NaN or Inf values.');
    end

    diagnostics.canonical_autocorrelations = full_correlations;
    diagnostics.rejected_components = rejected;
    diagnostics.component_unmixing = eeg_basis * full_transform;
    diagnostics.component_mixing = pinv(diagnostics.component_unmixing);
    diagnostics.skipped_reason = "";

    if size(cleaned, 1) ~= num_samples || size(cleaned, 2) ~= num_channels
        error('get_eeg_bss_cca_artifact_removal:InvalidReconstructionSize', ...
            'BSS-CCA reconstruction changed the EEG dimensions.');
    end
end

function diagnostics = empty_window_diagnostics()
    diagnostics = struct( ...
        'canonical_autocorrelations', zeros(0, 1), ...
        'rejected_components', zeros(1, 0), ...
        'component_unmixing', zeros(0, 0), ...
        'component_mixing', zeros(0, 0), ...
        'effective_eeg_rank', 0, ...
        'reconstruction_error', 0, ...
        'skipped_reason', "");
end

function [basis, effective_rank] = right_singular_basis(matrix)
    [~, singular_values, right_vectors] = svd(matrix, 'econ');
    values = diag(singular_values);
    if isempty(values) || max(values) <= eps
        basis = zeros(size(matrix, 2), 0);
        effective_rank = 0;
        return;
    end
    tolerance = max(size(matrix)) * eps(max(values));
    effective_rank = sum(values > tolerance);
    basis = right_vectors(:, 1:effective_rank);
end

function weight = overlap_weight(num_samples)
    phase = pi * ((1:num_samples).' - 0.5) / num_samples;
    weight = sin(phase).^2;
end

function metrics = calculate_metrics( ...
        original, cleaned, artifact, channel_names, delay_samples)
    num_channels = numel(channel_names);
    original_rms = zeros(num_channels, 1);
    cleaned_rms = zeros(num_channels, 1);
    artifact_rms = zeros(num_channels, 1);
    artifact_ratio = zeros(num_channels, 1);
    preservation_correlation = zeros(num_channels, 1);
    autocorrelation_before = zeros(num_channels, 1);
    autocorrelation_after = zeros(num_channels, 1);

    for channel_index = 1:num_channels
        before = center_signal(original(:, channel_index));
        after = center_signal(cleaned(:, channel_index));
        removed = artifact(:, channel_index);
        original_rms(channel_index) = root_mean_square(before);
        cleaned_rms(channel_index) = root_mean_square(after);
        artifact_rms(channel_index) = root_mean_square(removed);
        artifact_ratio(channel_index) = safe_ratio( ...
            artifact_rms(channel_index), original_rms(channel_index));
        preservation_correlation(channel_index) = ...
            absolute_correlation(before, after);
        autocorrelation_before(channel_index) = lagged_autocorrelation( ...
            before, delay_samples);
        autocorrelation_after(channel_index) = lagged_autocorrelation( ...
            after, delay_samples);
    end

    rms_change_percent = 100 * arrayfun(@(after, before) ...
        safe_ratio(after - before, before), cleaned_rms, original_rms);
    metrics = table(channel_names, original_rms, cleaned_rms, artifact_rms, ...
        artifact_ratio, preservation_correlation, autocorrelation_before, ...
        autocorrelation_after, rms_change_percent, ...
        'VariableNames', {'Channel', 'OriginalRMS', 'CleanedRMS', ...
        'ArtifactRMS', 'ArtifactToSignalRMSRatio', ...
        'OriginalToCleanedCorrelation', 'LagAutocorrelationBefore', ...
        'LagAutocorrelationAfter', 'RMSChangePercent'});
end

function value = lagged_autocorrelation(signal, delay_samples)
    if numel(signal) <= delay_samples
        value = 0;
    else
        value = absolute_correlation(signal(delay_samples + 1:end), ...
            signal(1:end - delay_samples));
    end
end

function centered = center_signal(signal)
    centered = double(signal(:));
    centered = centered - mean(centered);
end

function value = root_mean_square(signal)
    value = sqrt(mean(double(signal(:)).^2));
end

function value = safe_ratio(numerator, denominator)
    if denominator <= eps(max(1, abs(numerator)))
        if abs(numerator) <= eps
            value = 0;
        else
            value = Inf;
        end
    else
        value = numerator / denominator;
    end
end

function value = absolute_correlation(first_signal, second_signal)
    first_signal = center_signal(first_signal);
    second_signal = center_signal(second_signal);
    denominator = norm(first_signal) * norm(second_signal);
    if denominator <= eps
        value = 0;
    else
        value = abs((first_signal' * second_signal) / denominator);
    end
end
