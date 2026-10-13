function filtered = apply_eeg_filter_chain(data, Fs, settings)
%APPLY_EEG_FILTER_CHAIN Apply optional zero-phase band-pass and notch filters.
%   FILTERED = APPLY_EEG_FILTER_CHAIN(DATA, FS, SETTINGS) accepts either a
%   numeric vector/matrix (samples by channels) or a scalar structure of
%   channel vectors. SETTINGS fields are bandpass_enabled, band_hz,
%   bandpass_order, and notch_hz (0, 50, or 60).

    validateattributes(Fs, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'Fs', 2);
    settings = normalize_settings(settings);
    [matrix, converter] = to_matrix(data);

    if settings.bandpass_enabled
        band = double(settings.band_hz(:).');
        validateattributes(band, {'numeric'}, ...
            {'real','finite','vector','numel',2,'positive'}, mfilename, 'band_hz');
        if band(1) >= band(2) || band(2) >= Fs / 2
            error('apply_eeg_filter_chain:InvalidBand', ...
                'Band-pass edges must satisfy 0 < low < high < Nyquist (%.6g Hz).', Fs/2);
        end
        validateattributes(settings.bandpass_order, {'numeric'}, ...
            {'real','finite','scalar','integer','positive'}, mfilename, 'bandpass_order');
        if size(matrix, 1) <= 3 * settings.bandpass_order
            error('apply_eeg_filter_chain:SegmentTooShortForBandpass', ...
                'The segment is too short for zero-phase FIR order %d filtering.', ...
                settings.bandpass_order);
        end
        filterObject = designfilt('bandpassfir', ...
            'FilterOrder', settings.bandpass_order, ...
            'CutoffFrequency1', band(1), ...
            'CutoffFrequency2', band(2), ...
            'SampleRate', Fs);
        matrix = filtfilt(filterObject, matrix);
    end

    if settings.notch_hz ~= 0
        center = settings.notch_hz;
        if ~ismember(center, [50 60])
            error('apply_eeg_filter_chain:InvalidNotch', ...
                'notch_hz must be 0, 50, or 60 Hz.');
        end
        edges = [center - 1, center + 1];
        if edges(2) >= Fs / 2
            error('apply_eeg_filter_chain:NotchExceedsNyquist', ...
                'The %.0f Hz notch requires Nyquist above %.0f Hz.', center, edges(2));
        end
        [b, a] = butter(2, edges / (Fs / 2), 'stop');
        minimumLength = 3 * (max(numel(a), numel(b)) - 1) + 1;
        if size(matrix, 1) < minimumLength
            error('apply_eeg_filter_chain:SegmentTooShortForNotch', ...
                'The segment is too short for zero-phase notch filtering.');
        end
        matrix = filtfilt(b, a, matrix);
    end

    filtered = converter(matrix);
end

function settings = normalize_settings(settings)
    if nargin == 0 || isempty(settings)
        settings = struct();
    end
    if ~isstruct(settings) || ~isscalar(settings)
        error('apply_eeg_filter_chain:InvalidSettings', ...
            'settings must be a scalar structure.');
    end
    defaults = struct('bandpass_enabled', false, 'band_hz', [25 50], ...
        'bandpass_order', 20, 'notch_hz', 0);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(settings, names{index}) || isempty(settings.(names{index}))
            settings.(names{index}) = defaults.(names{index});
        end
    end
    settings.bandpass_enabled = logical(settings.bandpass_enabled);
    validateattributes(settings.notch_hz, {'numeric'}, ...
        {'real','finite','scalar'}, mfilename, 'notch_hz');
end

function [matrix, converter] = to_matrix(data)
    if isnumeric(data)
        if ~isreal(data) || isempty(data) || any(~isfinite(data(:)))
            error('apply_eeg_filter_chain:InvalidData', ...
                'Numeric data must be finite, real, and nonempty.');
        end
        wasRow = isrow(data);
        matrix = double(data);
        if wasRow
            matrix = matrix(:);
        end
        converter = @(x) restore_numeric(x, wasRow);
        return;
    end

    if ~isstruct(data) || ~isscalar(data)
        error('apply_eeg_filter_chain:InvalidData', ...
            'data must be numeric or a scalar channel structure.');
    end
    names = fieldnames(data);
    if isempty(names)
        error('apply_eeg_filter_chain:InvalidData', 'Channel structure is empty.');
    end
    numSamples = numel(data.(names{1}));
    matrix = zeros(numSamples, numel(names));
    orientations = false(numel(names), 1);
    for index = 1:numel(names)
        signal = data.(names{index});
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                numel(signal) ~= numSamples || any(~isfinite(signal(:)))
            error('apply_eeg_filter_chain:InvalidData', ...
                'All channel fields must be equal-length finite real vectors.');
        end
        orientations(index) = isrow(signal);
        matrix(:, index) = double(signal(:));
    end
    converter = @(x) restore_struct(x, names, orientations);
end

function value = restore_numeric(matrix, wasRow)
    value = matrix;
    if wasRow
        value = value.';
    end
end

function value = restore_struct(matrix, names, orientations)
    value = struct();
    for index = 1:numel(names)
        signal = matrix(:, index);
        if orientations(index)
            signal = signal.';
        end
        value.(names{index}) = signal;
    end
end
