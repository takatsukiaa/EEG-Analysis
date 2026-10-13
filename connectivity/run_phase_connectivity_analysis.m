function results = run_phase_connectivity_analysis(recording, options)
%RUN_PHASE_CONNECTIVITY_ANALYSIS Run aligned PLV, PLI, and wPLI analysis.
%   RESULTS = RUN_PHASE_CONNECTIVITY_ANALYSIS(RECORDING, OPTIONS) extracts
%   the selected half-open interval [start_time_sec,end_time_sec), applies
%   the configured in-memory filter chain, and calculates whole-segment and
%   sliding-window connectivity with GET_PHASE_CONNECTIVITY.

    validate_recording(recording);
    options = normalize_options(recording, options);
    selectedNames = resolve_channel_names(recording.channel_names, ...
        options.channel_names);
    if numel(selectedNames) < 2
        error('run_phase_connectivity_analysis:NotEnoughChannels', ...
            'Select at least two EEG channels for connectivity analysis.');
    end

    [startSample, endSample, actualStart, actualEnd] = interval_samples( ...
        recording, options.start_time_sec, options.end_time_sec);
    segment = struct();
    for index = 1:numel(selectedNames)
        name = selectedNames{index};
        signal = recording.data.(name);
        segment.(name) = double(signal(startSample:endSample));
    end

    filtered = apply_eeg_filter_chain(segment, recording.Fs, ...
        options.filter_settings);
    connectivity = get_phase_connectivity(filtered, recording.Fs, ...
        options.window_length_sec, options.window_step_sec);

    results = connectivity;
    results.absolute_time_axis = actualStart + connectivity.time_axis;
    results.start_time_sec = actualStart;
    results.end_time_sec = actualEnd;
    results.start_sample = startSample;
    results.end_sample = endSample;
    results.Fs = recording.Fs;
    results.source_metadata = recording.source_metadata;
    results.filter_settings = options.filter_settings;
    results.window_length_sec = options.window_length_sec;
    results.window_step_sec = options.window_step_sec;
    results.estimated_sliding_bytes = 3 * numel(connectivity.plv_3d_matrix) * 8;
end

function options = normalize_options(recording, options)
    if nargin < 2 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('run_phase_connectivity_analysis:InvalidOptions', ...
            'options must be a scalar structure.');
    end
    defaults = struct( ...
        'start_time_sec', 0, ...
        'end_time_sec', min(10, recording.duration_sec), ...
        'channel_names', {recording.channel_names}, ...
        'window_length_sec', 0.5, ...
        'window_step_sec', 0.1, ...
        'filter_settings', struct('bandpass_enabled', true, ...
            'band_hz', [25 50], 'bandpass_order', 20, 'notch_hz', 0));
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index}) || isempty(options.(names{index}))
            options.(names{index}) = defaults.(names{index});
        end
    end
    validateattributes(options.window_length_sec, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'window_length_sec');
    validateattributes(options.window_step_sec, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'window_step_sec');
end

function names = resolve_channel_names(allNames, requested)
    requested = cellstr(string(requested));
    names = cell(size(requested));
    for index = 1:numel(requested)
        match = find(strcmp(allNames, requested{index}), 1);
        if isempty(match)
            match = find(strcmpi(allNames, requested{index}), 1);
        end
        if isempty(match)
            error('run_phase_connectivity_analysis:ChannelNotFound', ...
                'Channel "%s" was not found in the recording.', requested{index});
        end
        names{index} = allNames{match};
    end
    if numel(unique(lower(string(names)))) ~= numel(names)
        error('run_phase_connectivity_analysis:DuplicateChannels', ...
            'Selected channel names must be unique.');
    end
end

function [firstSample, lastSample, actualStart, actualEnd] = ...
        interval_samples(recording, startTime, endTime)
    validateattributes(startTime, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'}, mfilename, 'start_time_sec');
    validateattributes(endTime, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'end_time_sec');
    if endTime <= startTime
        error('run_phase_connectivity_analysis:InvalidInterval', ...
            'End time must be greater than start time.');
    end
    if endTime > recording.duration_sec + eps(recording.duration_sec)
        error('run_phase_connectivity_analysis:IntervalExceedsRecording', ...
            'End time %.6g exceeds recording duration %.6g seconds.', ...
            endTime, recording.duration_sec);
    end
    firstSample = floor(startTime * recording.Fs) + 1;
    lastSample = min(floor(endTime * recording.Fs), recording.num_samples);
    if lastSample - firstSample + 1 < 2
        error('run_phase_connectivity_analysis:IntervalTooShort', ...
            'The selected interval contains fewer than two samples.');
    end
    actualStart = (firstSample - 1) / recording.Fs;
    actualEnd = lastSample / recording.Fs;
end

function validate_recording(recording)
    required = {'data','Fs','channel_names','num_samples','duration_sec','source_metadata'};
    if ~isstruct(recording) || ~isscalar(recording) || ...
            ~all(isfield(recording, required))
        error('run_phase_connectivity_analysis:InvalidRecording', ...
            'recording must be returned by load_phase_connectivity_recording.');
    end
end
