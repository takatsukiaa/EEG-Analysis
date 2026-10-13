function result = compute_single_channel_psd(recording, options)
%COMPUTE_SINGLE_CHANNEL_PSD Calculate raw and independently filtered PSD.

    required = {'data','Fs','channel_names','num_samples','duration_sec','source_metadata'};
    if ~isstruct(recording) || ~isscalar(recording) || ...
            ~all(isfield(recording, required))
        error('compute_single_channel_psd:InvalidRecording', ...
            'recording must be returned by load_phase_connectivity_recording.');
    end
    options = normalize_options(recording, options);
    channelName = resolve_channel(recording.channel_names, options.channel_name);

    [firstSample, lastSample] = interval_samples(recording, ...
        options.start_time_sec, options.end_time_sec);
    rawSignal = double(recording.data.(channelName)(firstSample:lastSample));
    filteredSignal = apply_eeg_filter_chain(rawSignal, recording.Fs, ...
        options.filter_settings);

    rawSignal = rawSignal(:) - mean(rawSignal, 'omitnan');
    filteredSignal = filteredSignal(:) - mean(filteredSignal, 'omitnan');
    windowSamples = max(4, round(options.window_length_sec * recording.Fs));
    windowSamples = min(windowSamples, numel(rawSignal));
    overlapSamples = floor(windowSamples * options.overlap_fraction);
    nfft = max(256, 2 ^ nextpow2(windowSamples));
    window = hamming(windowSamples, 'periodic');

    [rawPsd, frequency] = pwelch(rawSignal, window, overlapSamples, ...
        nfft, recording.Fs);
    filteredPsd = pwelch(filteredSignal, window, overlapSamples, ...
        nfft, recording.Fs);

    maximumFrequency = min(options.max_frequency_hz, recording.Fs / 2);
    keep = frequency <= maximumFrequency;
    result = struct();
    result.channel_name = channelName;
    result.frequency = frequency(keep);
    result.raw_psd = rawPsd(keep);
    result.filtered_psd = filteredPsd(keep);
    result.Fs = recording.Fs;
    result.start_sample = firstSample;
    result.end_sample = lastSample;
    result.start_time_sec = (firstSample - 1) / recording.Fs;
    result.end_time_sec = lastSample / recording.Fs;
    result.window_length_sec = windowSamples / recording.Fs;
    result.overlap_fraction = options.overlap_fraction;
    result.nfft = nfft;
    result.max_frequency_hz = maximumFrequency;
    result.filter_settings = options.filter_settings;
    result.source_metadata = recording.source_metadata;
end

function options = normalize_options(recording, options)
    if nargin < 2 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('compute_single_channel_psd:InvalidOptions', ...
            'options must be a scalar structure.');
    end
    defaults = struct( ...
        'start_time_sec', 0, ...
        'end_time_sec', min(10, recording.duration_sec), ...
        'channel_name', recording.channel_names{1}, ...
        'window_length_sec', 1, ...
        'overlap_fraction', 0.5, ...
        'max_frequency_hz', 70, ...
        'filter_settings', struct('bandpass_enabled', false, ...
            'band_hz', [1 70], 'bandpass_order', 20, 'notch_hz', 0));
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index}) || isempty(options.(names{index}))
            options.(names{index}) = defaults.(names{index});
        end
    end
    validateattributes(options.window_length_sec, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'window_length_sec');
    validateattributes(options.overlap_fraction, {'numeric'}, ...
        {'real','finite','scalar','>=',0,'<',1}, mfilename, 'overlap_fraction');
    validateattributes(options.max_frequency_hz, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'max_frequency_hz');
end

function name = resolve_channel(allNames, requested)
    requested = char(string(requested));
    match = find(strcmp(allNames, requested), 1);
    if isempty(match)
        match = find(strcmpi(allNames, requested), 1);
    end
    if isempty(match)
        error('compute_single_channel_psd:ChannelNotFound', ...
            'Channel "%s" was not found in the recording.', requested);
    end
    name = allNames{match};
end

function [firstSample, lastSample] = interval_samples(recording, startTime, endTime)
    validateattributes(startTime, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'}, mfilename, 'start_time_sec');
    validateattributes(endTime, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'end_time_sec');
    if endTime <= startTime || endTime > recording.duration_sec + eps(recording.duration_sec)
        error('compute_single_channel_psd:InvalidInterval', ...
            'Use an interval satisfying 0 <= start < end <= %.6g seconds.', ...
            recording.duration_sec);
    end
    firstSample = floor(startTime * recording.Fs) + 1;
    lastSample = min(floor(endTime * recording.Fs), recording.num_samples);
    if lastSample - firstSample + 1 < 4
        error('compute_single_channel_psd:IntervalTooShort', ...
            'The selected interval contains fewer than four samples.');
    end
end
