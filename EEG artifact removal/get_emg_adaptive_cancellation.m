function results = get_emg_adaptive_cancellation(eeg_data, emg_reference, Fs, options)
%GET_EMG_ADAPTIVE_CANCELLATION Remove EMG-correlated activity from EEG.
%   RESULTS = GET_EMG_ADAPTIVE_CANCELLATION(EEG_DATA, EMG_REFERENCE, FS)
%   applies one normalized-LMS adaptive FIR filter per EEG channel. EEG_DATA
%   is a scalar structure whose fields are channel vectors. EMG_REFERENCE is
%   a synchronized reference vector with the same number of samples.
%
%   Name-value options:
%       FilterLengthSec - adaptive FIR duration in seconds (default 0.08)
%       StepSize        - normalized-LMS step size in (0, 2) (default 0.05)
%
%   The adaptive-filter output estimates the EMG-related artifact and the
%   error output is returned as cleaned EEG. Diagnostic metrics quantify
%   signal changes; they do not establish physiological validity.

    arguments
        eeg_data (1, 1) struct
        emg_reference {mustBeNumeric, mustBeReal, mustBeVector}
        Fs (1, 1) double {mustBeFinite, mustBePositive}
        options.FilterLengthSec (1, 1) double {mustBeFinite, mustBePositive} = 0.08
        options.StepSize (1, 1) double {mustBeFinite, mustBePositive} = 0.05
    end

    if options.StepSize >= 2
        error('get_emg_adaptive_cancellation:InvalidStepSize', ...
            'StepSize must be greater than zero and less than 2 for normalized LMS.');
    end
    if exist('dsp.LMSFilter', 'class') ~= 8
        error('get_emg_adaptive_cancellation:MissingDSPToolbox', ...
            ['DSP System Toolbox is required. The class dsp.LMSFilter ' ...
             'was not found on the MATLAB path.']);
    end

    channel_names = fieldnames(eeg_data);
    if isempty(channel_names)
        error('get_emg_adaptive_cancellation:NoEEGChannels', ...
            'eeg_data must contain at least one EEG channel.');
    end

    emg_reference = double(emg_reference(:));
    if isempty(emg_reference) || any(~isfinite(emg_reference))
        error('get_emg_adaptive_cancellation:InvalidEMGReference', ...
            'emg_reference must be a nonempty finite numeric vector.');
    end
    num_samples = numel(emg_reference);

    eeg_matrix = zeros(num_samples, numel(channel_names));
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = eeg_data.(channel_name);
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('get_emg_adaptive_cancellation:InvalidEEGChannel', ...
                'EEG channel "%s" must be a nonempty finite real numeric vector.', ...
                channel_name);
        end
        if numel(signal) ~= num_samples
            error('get_emg_adaptive_cancellation:UnequalChannelLengths', ...
                ['EEG channel "%s" has %d samples, but the EMG reference ' ...
                 'has %d samples.'], channel_name, numel(signal), num_samples);
        end
        eeg_matrix(:, channel_index) = double(signal(:));
    end

    filter_length_samples = max(1, round(options.FilterLengthSec * Fs));
    if num_samples <= filter_length_samples
        error('get_emg_adaptive_cancellation:SignalTooShort', ...
            ['The signal has %d samples; it must contain more samples than ' ...
             'the %d-tap adaptive filter.'], num_samples, filter_length_samples);
    end

    % A median baseline preserves exact quiet/zero regions when the reference
    % contains sparse bursts. Mean centering would turn those regions into a
    % constant nonzero filter input and could create an artificial EEG shift.
    emg_baseline = median(emg_reference);
    normalized_emg = emg_reference - emg_baseline;
    emg_rms = root_mean_square(normalized_emg);
    if emg_rms <= eps(max(1, max(abs(emg_reference))))
        error('get_emg_adaptive_cancellation:ConstantEMGReference', ...
            'The EMG reference has no usable variation after demeaning.');
    end
    activity_threshold = max(eps(max(1, max(abs(normalized_emg)))), ...
        1e-6 * max(abs(normalized_emg)));
    emg_activity_fraction = mean(abs(normalized_emg) > activity_threshold);
    if emg_activity_fraction < 0.01
        warning('get_emg_adaptive_cancellation:SparseEMGReference', ...
            ['Only %.3g%% of EMG-reference samples differ from its median. ' ...
             'Verify that this channel contains a usable EMG waveform.'], ...
            100 * emg_activity_fraction);
    end
    normalized_emg = normalized_emg / emg_rms;

    eeg_means = mean(eeg_matrix, 1);
    demeaned_eeg = eeg_matrix - eeg_means;
    cleaned_matrix = zeros(size(eeg_matrix));
    artifact_matrix = zeros(size(eeg_matrix));
    final_weights = zeros(filter_length_samples, numel(channel_names));

    for channel_index = 1:numel(channel_names)
        adaptive_filter = dsp.LMSFilter(filter_length_samples, ...
            'Method', 'Normalized LMS', ...
            'StepSize', options.StepSize);
        [artifact, residual, weights] = adaptive_filter( ...
            normalized_emg, demeaned_eeg(:, channel_index));

        artifact_matrix(:, channel_index) = double(artifact);
        cleaned_matrix(:, channel_index) = double(residual) + eeg_means(channel_index);
        final_weights(:, channel_index) = double(weights(:));
    end

    cleaned_eeg = struct();
    artifact_estimate = struct();
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        cleaned_eeg.(channel_name) = cleaned_matrix(:, channel_index);
        artifact_estimate.(channel_name) = artifact_matrix(:, channel_index);
    end

    original_rms = zeros(numel(channel_names), 1);
    cleaned_rms = zeros(numel(channel_names), 1);
    artifact_rms = zeros(numel(channel_names), 1);
    artifact_to_signal_ratio = zeros(numel(channel_names), 1);
    emg_correlation_before = zeros(numel(channel_names), 1);
    emg_correlation_after = zeros(numel(channel_names), 1);

    for channel_index = 1:numel(channel_names)
        original_centered = eeg_matrix(:, channel_index) - eeg_means(channel_index);
        cleaned_centered = cleaned_matrix(:, channel_index) - ...
            mean(cleaned_matrix(:, channel_index));
        original_rms(channel_index) = root_mean_square(original_centered);
        cleaned_rms(channel_index) = root_mean_square(cleaned_centered);
        artifact_rms(channel_index) = root_mean_square(artifact_matrix(:, channel_index));
        artifact_to_signal_ratio(channel_index) = safe_ratio( ...
            artifact_rms(channel_index), original_rms(channel_index));
        emg_correlation_before(channel_index) = absolute_correlation( ...
            original_centered, normalized_emg);
        emg_correlation_after(channel_index) = absolute_correlation( ...
            cleaned_centered, normalized_emg);
    end

    rms_change_percent = 100 * arrayfun(@(cleaned, original) ...
        safe_ratio(cleaned - original, original), cleaned_rms, original_rms);

    metrics = table(channel_names, original_rms, cleaned_rms, artifact_rms, ...
        artifact_to_signal_ratio, emg_correlation_before, ...
        emg_correlation_after, rms_change_percent, ...
        'VariableNames', {'Channel', 'OriginalRMS', 'CleanedRMS', ...
        'ArtifactRMS', 'ArtifactToSignalRMSRatio', ...
        'EMGCorrelationBefore', 'EMGCorrelationAfter', 'RMSChangePercent'});

    results = struct();
    results.cleaned_eeg = cleaned_eeg;
    results.artifact_estimate = artifact_estimate;
    results.channel_names = channel_names;
    results.final_weights = final_weights;
    results.filter_length_sec = options.FilterLengthSec;
    results.filter_length_samples = filter_length_samples;
    results.step_size = options.StepSize;
    results.emg_reference_baseline = emg_baseline;
    results.emg_reference_rms = emg_rms;
    results.emg_reference_activity_fraction = emg_activity_fraction;
    results.metrics = metrics;
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
    first_signal = double(first_signal(:));
    second_signal = double(second_signal(:));
    first_signal = first_signal - mean(first_signal);
    second_signal = second_signal - mean(second_signal);
    denominator = norm(first_signal) * norm(second_signal);
    if denominator <= eps
        value = 0;
    else
        value = abs((first_signal' * second_signal) / denominator);
    end
end
