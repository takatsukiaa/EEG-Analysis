function result = compute_psd_topography(recording, options)
%COMPUTE_PSD_TOPOGRAPHY Calculate multichannel band power by epoch.
%   RESULT = COMPUTE_PSD_TOPOGRAPHY(RECORDING, OPTIONS) analyzes every EEG
%   channel in consecutive, non-overlapping epochs. Raw and independently
%   filtered band powers are retained; waveform data are never returned.

    required = {'data','Fs','channel_names','num_samples','duration_sec','source_metadata'};
    if ~isstruct(recording) || ~isscalar(recording) || ...
            ~all(isfield(recording, required))
        error('compute_psd_topography:InvalidRecording', ...
            'recording must be returned by load_phase_connectivity_recording.');
    end
    options = normalize_options(recording, options);
    validate_band(options.band_hz, recording.Fs);

    [requestedFirst, requestedLast] = interval_samples(recording, ...
        options.start_time_sec, options.end_time_sec);
    epochSamples = round(options.epoch_length_sec * recording.Fs);
    if epochSamples < 4
        error('compute_psd_topography:EpochTooShort', ...
            'Each epoch must contain at least four samples.');
    end
    requestedSamples = requestedLast - requestedFirst + 1;
    numberOfEpochs = floor(requestedSamples / epochSamples);
    if numberOfEpochs < 1
        error('compute_psd_topography:IntervalTooShort', ...
            'The selected interval must contain at least one complete %.6g-second epoch.', ...
            epochSamples / recording.Fs);
    end

    analyzedSamples = numberOfEpochs * epochSamples;
    analyzedLast = requestedFirst + analyzedSamples - 1;
    unusedSamples = requestedSamples - analyzedSamples;
    channelNames = cellstr(string(recording.channel_names));
    segment = zeros(analyzedSamples, numel(channelNames));
    for channelIndex = 1:numel(channelNames)
        signal = recording.data.(channelNames{channelIndex});
        segment(:, channelIndex) = double(signal(requestedFirst:analyzedLast));
    end
    filteredSegment = apply_eeg_filter_chain(segment, recording.Fs, ...
        options.filter_settings);

    nfft = max(256, 2 ^ nextpow2(epochSamples));
    epochWindow = hamming(epochSamples, 'periodic');
    rawBandPower = zeros(numberOfEpochs, numel(channelNames));
    filteredBandPower = zeros(numberOfEpochs, numel(channelNames));
    epochStartTime = zeros(numberOfEpochs, 1);
    epochEndTime = zeros(numberOfEpochs, 1);
    epochCenterTime = zeros(numberOfEpochs, 1);

    for epochIndex = 1:numberOfEpochs
        first = (epochIndex - 1) * epochSamples + 1;
        last = first + epochSamples - 1;
        rawEpoch = detrend(segment(first:last, :), 0);
        filteredEpoch = detrend(filteredSegment(first:last, :), 0);
        [rawPsd, frequency] = pwelch(rawEpoch, epochWindow, 0, nfft, recording.Fs);
        filteredPsd = pwelch(filteredEpoch, epochWindow, 0, nfft, recording.Fs);
        rawBandPower(epochIndex, :) = bandpower( ...
            rawPsd, frequency, options.band_hz, 'psd');
        filteredBandPower(epochIndex, :) = bandpower( ...
            filteredPsd, frequency, options.band_hz, 'psd');

        absoluteFirst = requestedFirst + first - 1;
        absoluteLast = requestedFirst + last - 1;
        epochStartTime(epochIndex) = (absoluteFirst - 1) / recording.Fs;
        epochEndTime(epochIndex) = absoluteLast / recording.Fs;
        epochCenterTime(epochIndex) = ...
            (epochStartTime(epochIndex) + epochEndTime(epochIndex)) / 2;
    end

    result = struct();
    result.channel_names = channelNames;
    result.Fs = recording.Fs;
    result.band_hz = double(options.band_hz(:).');
    result.filter_settings = options.filter_settings;
    result.requested_start_sample = requestedFirst;
    result.requested_end_sample = requestedLast;
    result.start_sample = requestedFirst;
    result.end_sample = analyzedLast;
    result.start_time_sec = (requestedFirst - 1) / recording.Fs;
    result.end_time_sec = analyzedLast / recording.Fs;
    result.requested_end_time_sec = requestedLast / recording.Fs;
    result.epoch_length_sec = epochSamples / recording.Fs;
    result.epoch_samples = epochSamples;
    result.num_epochs = numberOfEpochs;
    result.epoch_start_time_sec = epochStartTime;
    result.epoch_end_time_sec = epochEndTime;
    result.epoch_center_time_sec = epochCenterTime;
    result.unused_tail_samples = unusedSamples;
    result.unused_tail_sec = unusedSamples / recording.Fs;
    result.nfft = nfft;
    result.raw_band_power = rawBandPower;
    result.filtered_band_power = filteredBandPower;
    result.raw_whole_band_power = mean(rawBandPower, 1);
    result.filtered_whole_band_power = mean(filteredBandPower, 1);
    result.source_metadata = recording.source_metadata;
end

function options = normalize_options(recording, options)
    if nargin < 2 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('compute_psd_topography:InvalidOptions', ...
            'options must be a scalar structure.');
    end
    defaults = struct( ...
        'start_time_sec', 0, ...
        'end_time_sec', min(10, recording.duration_sec), ...
        'band_hz', [25 50], ...
        'epoch_length_sec', 0.5, ...
        'filter_settings', struct('bandpass_enabled', false, ...
            'band_hz', [1 70], 'bandpass_order', 20, 'notch_hz', 0));
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index}) || isempty(options.(names{index}))
            options.(names{index}) = defaults.(names{index});
        end
    end
    validateattributes(options.epoch_length_sec, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'epoch_length_sec');
end

function validate_band(band, Fs)
    validateattributes(band, {'numeric'}, ...
        {'real','finite','vector','numel',2,'nonnegative'}, mfilename, 'band_hz');
    band = double(band(:).');
    if band(1) >= band(2) || band(2) > Fs / 2
        error('compute_psd_topography:InvalidBand', ...
            'Band edges must satisfy 0 <= low < high <= Nyquist (%.6g Hz).', Fs/2);
    end
end

function [firstSample, lastSample] = interval_samples(recording, startTime, endTime)
    validateattributes(startTime, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'}, mfilename, 'start_time_sec');
    validateattributes(endTime, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'end_time_sec');
    if endTime <= startTime || endTime > recording.duration_sec + eps(recording.duration_sec)
        error('compute_psd_topography:InvalidInterval', ...
            'Use an interval satisfying 0 <= start < end <= %.6g seconds.', ...
            recording.duration_sec);
    end
    firstSample = floor(startTime * recording.Fs) + 1;
    lastSample = min(floor(endTime * recording.Fs), recording.num_samples);
end
