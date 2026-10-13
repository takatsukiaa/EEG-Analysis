function recording = load_phase_connectivity_recording(inputFile)
%LOAD_PHASE_CONNECTIVITY_RECORDING Load an EDF or MAT EEG recording.
%   RECORDING = LOAD_PHASE_CONNECTIVITY_RECORDING(INPUTFILE) returns a
%   normalized recording structure with scalar fields Fs, num_samples and
%   duration_sec, a channel-vector structure in data, channel_names, source
%   metadata, and optional EEGLAB-compatible chanlocs.

    inputFile = char(string(inputFile));
    if isempty(inputFile) || ~isfile(inputFile)
        error('load_phase_connectivity_recording:FileNotFound', ...
            'EEG file was not found: %s', inputFile);
    end

    [~, fileStem, extension] = fileparts(inputFile);
    switch lower(extension)
        case '.edf'
            [data, Fs, chanlocs] = load_edf(inputFile);
        case '.mat'
            [data, Fs, chanlocs] = load_mat(inputFile);
        otherwise
            error('load_phase_connectivity_recording:UnsupportedFileType', ...
                'Unsupported extension "%s". Select an EDF or MAT file.', extension);
    end

    data = order_channels(data);
    [numSamples, channelNames] = validate_channels(data);

    recording = struct();
    recording.input_file = inputFile;
    recording.file_name = [fileStem extension];
    recording.extension = lower(extension);
    recording.data = data;
    recording.Fs = Fs;
    recording.channel_names = channelNames;
    recording.num_samples = numSamples;
    recording.duration_sec = numSamples / Fs;
    recording.chanlocs = chanlocs;
    recording.source_metadata = struct( ...
        'input_file', inputFile, ...
        'file_name', [fileStem extension], ...
        'extension', lower(extension), ...
        'sampling_rate_hz', Fs, ...
        'num_samples', numSamples, ...
        'duration_sec', numSamples / Fs, ...
        'num_channels', numel(channelNames));
end

function [data, Fs, chanlocs] = load_edf(inputFile)
    info = edfinfo(inputFile);
    recordDuration = seconds(info.DataRecordDuration);
    if ~isscalar(recordDuration) || ~isfinite(recordDuration) || recordDuration <= 0
        error('load_phase_connectivity_recording:InvalidEDFRecordDuration', ...
            'EDF data-record duration must be a positive scalar.');
    end

    sampleRates = double(info.NumSamples(:)) ./ recordDuration;
    tableData = edfread(inputFile);
    variableNames = tableData.Properties.VariableNames;
    data = struct();
    retainedRates = [];

    for index = 1:numel(variableNames)
        name = variableNames{index};
        if is_non_eeg_name(name) || index > numel(sampleRates)
            continue;
        end
        value = tableData.(name);
        if iscell(value)
            if isempty(value) || any(~cellfun(@isnumeric, value))
                continue;
            end
            blocks = cellfun(@(x) double(x(:)), value, 'UniformOutput', false);
            signal = vertcat(blocks{:});
        elseif isnumeric(value)
            signal = double(value(:));
        else
            continue;
        end
        validate_signal(signal, name);
        data.(name) = signal;
        retainedRates(end + 1) = sampleRates(index); %#ok<AGROW>
    end

    if isempty(retainedRates)
        error('load_phase_connectivity_recording:NoEEGChannels', ...
            'No numeric EEG channels were found in the EDF file.');
    end
    Fs = retainedRates(1);
    tolerance = max(1, abs(Fs)) * 1e-10;
    if any(abs(retainedRates - Fs) > tolerance)
        error('load_phase_connectivity_recording:MixedSamplingRates', ...
            'Retained EDF channels have different sampling rates.');
    end
    chanlocs = struct([]);
end

function [data, Fs, chanlocs] = load_mat(inputFile)
    loaded = load(inputFile);
    chanlocs = extract_chanlocs(loaded);

    if isfield(loaded, 'EEG') && isstruct(loaded.EEG) && ...
            isfield(loaded.EEG, 'data')
        EEG = loaded.EEG;
        if isfield(EEG, 'srate')
            Fs = double(EEG.srate);
        else
            Fs = find_sample_rate(loaded);
        end
        [matrix, labels] = normalize_matrix(double(EEG.data), chanlocs);
        data = matrix_to_struct(matrix, labels);
        return;
    end

    Fs = find_sample_rate(loaded);
    ignored = {'Fs','fs','srate','sampleRate','samplingRate','SamplingRate', ...
        'data','EEG','chanlocs','urchanlocs'};
    names = fieldnames(loaded);
    vectorNames = {};
    vectorLengths = [];
    for index = 1:numel(names)
        name = names{index};
        value = loaded.(name);
        if any(strcmpi(name, ignored)) || is_non_eeg_name(name) || ...
                contains(name, '_segment', 'IgnoreCase', true)
            continue;
        end
        if isnumeric(value) && isreal(value) && isvector(value) && numel(value) > 1
            validate_signal(value, name);
            vectorNames{end + 1} = name; %#ok<AGROW>
            vectorLengths(end + 1) = numel(value); %#ok<AGROW>
        end
    end

    if ~isempty(vectorNames)
        targetLength = max(vectorLengths);
        keep = vectorLengths == targetLength;
        vectorNames = vectorNames(keep);
        data = struct();
        for index = 1:numel(vectorNames)
            name = vectorNames{index};
            data.(name) = double(loaded.(name)(:));
        end
        return;
    end

    if isfield(loaded, 'data') && isnumeric(loaded.data) && ismatrix(loaded.data)
        [matrix, labels] = normalize_matrix(double(loaded.data), chanlocs);
        data = matrix_to_struct(matrix, labels);
        return;
    end

    error('load_phase_connectivity_recording:NoEEGChannels', ...
        ['No EEG channels were detected. Expected named numeric vectors, ' ...
         'an EEG.data structure, or a numeric data matrix.']);
end

function Fs = find_sample_rate(loaded)
    candidates = {'Fs','fs','srate','sampleRate','samplingRate','SamplingRate'};
    for index = 1:numel(candidates)
        name = candidates{index};
        if isfield(loaded, name) && isnumeric(loaded.(name)) && ...
                isscalar(loaded.(name))
            Fs = double(loaded.(name));
            validateattributes(Fs, {'numeric'}, ...
                {'real','finite','scalar','positive'}, mfilename, name);
            return;
        end
    end
    error('load_phase_connectivity_recording:MissingSamplingRate', ...
        'MAT file must contain a scalar sampling rate such as Fs, fs, or srate.');
end

function chanlocs = extract_chanlocs(loaded)
    chanlocs = struct([]);
    if isfield(loaded, 'chanlocs') && isstruct(loaded.chanlocs)
        chanlocs = loaded.chanlocs;
    elseif isfield(loaded, 'EEG') && isstruct(loaded.EEG) && ...
            isfield(loaded.EEG, 'chanlocs') && isstruct(loaded.EEG.chanlocs)
        chanlocs = loaded.EEG.chanlocs;
    end
end

function [matrix, labels] = normalize_matrix(raw, chanlocs)
    if ndims(raw) > 2
        raw = reshape(raw, size(raw, 1), []);
    end
    if ~ismatrix(raw) || isempty(raw) || ~isreal(raw) || any(~isfinite(raw(:)))
        error('load_phase_connectivity_recording:InvalidDataMatrix', ...
            'EEG data matrix must be finite, real, and nonempty.');
    end

    numLocationLabels = numel(chanlocs);
    if numLocationLabels > 0 && size(raw, 1) == numLocationLabels
        matrix = raw.';
    elseif numLocationLabels > 0 && size(raw, 2) == numLocationLabels
        matrix = raw;
    elseif size(raw, 1) <= size(raw, 2)
        matrix = raw.';
    else
        matrix = raw;
    end

    numChannels = size(matrix, 2);
    if numLocationLabels == numChannels && isfield(chanlocs, 'labels')
        labels = cellstr(string({chanlocs.labels}));
    else
        labels = arrayfun(@(x) sprintf('Ch%d', x), 1:numChannels, ...
            'UniformOutput', false);
    end
end

function data = matrix_to_struct(matrix, labels)
    data = struct();
    labels = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(labels));
    for index = 1:size(matrix, 2)
        if is_non_eeg_name(labels{index})
            continue;
        end
        data.(labels{index}) = double(matrix(:, index));
    end
end

function data = order_channels(data)
    preferred = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','Cz', ...
        'C4','T4','T5','P3','Pz','P4','T6','O1','O2'};
    names = fieldnames(data);
    canonical = cellfun(@canonical_label, names, 'UniformOutput', false);
    used = false(size(names));
    ordered = {};
    for index = 1:numel(preferred)
        match = find(strcmpi(canonical, canonical_label(preferred{index})) & ~used, 1);
        if ~isempty(match)
            ordered{end + 1, 1} = names{match}; %#ok<AGROW>
            used(match) = true;
        end
    end
    ordered = [ordered; names(~used)];
    data = orderfields(data, ordered);
end

function [numSamples, channelNames] = validate_channels(data)
    channelNames = fieldnames(data);
    if isempty(channelNames)
        error('load_phase_connectivity_recording:NoEEGChannels', ...
            'No usable EEG channels were found.');
    end
    numSamples = numel(data.(channelNames{1}));
    for index = 1:numel(channelNames)
        signal = data.(channelNames{index});
        validate_signal(signal, channelNames{index});
        if numel(signal) ~= numSamples
            error('load_phase_connectivity_recording:UnequalChannelLengths', ...
                'Channel "%s" has %d samples; expected %d.', ...
                channelNames{index}, numel(signal), numSamples);
        end
    end
end

function validate_signal(signal, name)
    if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || isempty(signal) || ...
            any(~isfinite(signal(:)))
        error('load_phase_connectivity_recording:InvalidChannel', ...
            'Channel "%s" must be a finite, real numeric vector.', name);
    end
end

function tf = is_non_eeg_name(name)
    tf = contains(name, 'emg', 'IgnoreCase', true) || ...
         contains(name, 'ekg', 'IgnoreCase', true) || ...
         contains(name, 'ecg', 'IgnoreCase', true) || ...
         contains(name, 'annotation', 'IgnoreCase', true);
end

function value = canonical_label(value)
    value = upper(regexprep(char(string(value)), '[^A-Za-z0-9]', ''));
    value = regexprep(value, 'AVG$', '');
end
