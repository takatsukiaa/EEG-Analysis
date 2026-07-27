%% Plot one combined PSD figure for selected EEG channels
% Edit the settings in this section, then run the script.

clc;
clearvars -except filePath channelNames channelName channelIndex ...
    startTimeSec endTimeSec timePointSec secondsBefore secondsAfter ...
    frequencyRangeHz spectrogramWindowSec ...
    spectrogramOverlapFraction nfft minMaxNormalizeSpectrum ...
    showPowerDb timeTickFormat outputFolder outputResolutionDpi ...
    outputImageFormat showFigures closeFiguresAfterSave;
close all;

%% User settings
if ~exist('filePath', 'var')
    filePath = "";              % Example: "Raw_Data/No2.mat"; leave empty to choose a file
end
if ~exist('channelNames', 'var')
    channelNames = {"C4" "F8" "P4" "T4" "Pz"};
end
if ~exist('startTimeSec', 'var')
    startTimeSec = 26441;
end
if ~exist('endTimeSec', 'var')
    endTimeSec = 26446;
end
if ~exist('timePointSec', 'var')
    timePointSec = 26443.5;          % Example: 100; [] uses startTimeSec/endTimeSec instead
end
if ~exist('secondsBefore', 'var')
    secondsBefore = 2.5;
end
if ~exist('secondsAfter', 'var')
    secondsAfter = 2.5;
end

if ~exist('frequencyRangeHz', 'var')
    frequencyRangeHz = [0 70]; % Match the example-style view
end
if ~exist('windowLengthSec', 'var')
    windowLengthSec = 1;
end
if ~exist('overlapFraction', 'var')
    overlapFraction = 0.5;
end
if ~exist('plotInDb', 'var')
    plotInDb = false;           % false gives PSD in uV^2/Hz like the example
end
if ~exist('highlightFrequencyRangeHz', 'var')
    highlightFrequencyRangeHz = [25 50];   % [] disables the highlighted frequency band
end
if ~exist('highlightBandColor', 'var')
    highlightBandColor = [1.0 0.90 0.25];
end
if ~exist('highlightBandAlpha', 'var')
    highlightBandAlpha = 0.18;
end
if ~exist('highlightLineColor', 'var')
    highlightLineColor = [0.85 0.45 0.00];
end
if ~exist('outputFolder', 'var')
    outputFolder = fullfile(pwd, 'Yeh_PSDandSpectrogram/2大/26441/');
end
if ~exist('outputResolutionDpi', 'var')
    outputResolutionDpi = 200;
end
if ~exist('outputImageFormat', 'var')
    outputImageFormat = "png";
end
if ~exist('showFigures', 'var')
    showFigures = false;
end
if ~exist('closeFiguresAfterSave', 'var')
    closeFiguresAfterSave = true;
end

%% Select and load file
if strlength(string(filePath)) == 0
    [fileName, folderName] = uigetfile( ...
        {'*.edf;*.mat', 'EEG files (*.edf, *.mat)'; ...
         '*.edf', 'EDF files (*.edf)'; ...
         '*.mat', 'MAT files (*.mat)'}, ...
        'Select an EDF or MAT file');

    if isequal(fileName, 0)
        error('No file was selected.');
    end

    filePath = fullfile(folderName, fileName);
end

%% Determine requested time period
if isempty(timePointSec)
    requestedStartTimeSec = startTimeSec;
    requestedEndTimeSec = endTimeSec;
else
    if ~isscalar(timePointSec)
        error('timePointSec must be a scalar value in seconds.');
    end
    timePointSec = double(timePointSec);
    if timePointSec < 0
        error('timePointSec must be nonnegative.');
    end
    if ~isscalar(secondsBefore) || ~isscalar(secondsAfter)
        error('secondsBefore and secondsAfter must be scalar values in seconds.');
    end
    secondsBefore = double(secondsBefore);
    secondsAfter = double(secondsAfter);
    if secondsBefore < 0 || secondsAfter < 0 || secondsBefore + secondsAfter <= 0
        error('secondsBefore and secondsAfter must be nonnegative and cannot both be zero.');
    end

    requestedStartTimeSec = timePointSec - secondsBefore;
    requestedEndTimeSec = timePointSec + secondsAfter;
end

%% Load selected channels
[eegData, Fs, labels, dataStartTimeSec, recordDurationSec] = loadEegFileForPsd( ...
    filePath, channelNames, requestedStartTimeSec, requestedEndTimeSec);

%% Select time period
nSamples = size(eegData, 1);

if requestedStartTimeSec < 0
    warning('Requested start time %.3f sec is before the recording start. Using 0 sec.', ...
        requestedStartTimeSec);
end

if requestedEndTimeSec > recordDurationSec
    warning(['Requested end time %.3f sec is after the recording end %.3f sec. ', ...
        'Using the end of the recording.'], requestedEndTimeSec, recordDurationSec);
end

if requestedEndTimeSec <= 0 || requestedStartTimeSec >= recordDurationSec || ...
        requestedEndTimeSec <= requestedStartTimeSec
    error('Use a valid time period.');
end
if overlapFraction < 0 || overlapFraction >= 1
    error('overlapFraction must be >= 0 and < 1.');
end

firstSample = max(1, floor((requestedStartTimeSec - dataStartTimeSec) * Fs) + 1);
lastSample = min(nSamples, floor((requestedEndTimeSec - dataStartTimeSec) * Fs));
segment = eegData(firstSample:lastSample, :);
actualStartTimeSec = dataStartTimeSec + (firstSample - 1) / Fs;
actualEndTimeSec = dataStartTimeSec + lastSample / Fs;

if size(segment, 1) < 4
    error('Selected segment is too short for PSD calculation.');
end

%% Compute PSD for each channel
windowLengthSamples = max(4, round(windowLengthSec * Fs));
windowLengthSamples = min(windowLengthSamples, size(segment, 1));
overlapSamples = floor(windowLengthSamples * overlapFraction);
nfft = max(256, 2 ^ nextpow2(windowLengthSamples));
window = hamming(windowLengthSamples, 'periodic');

segmentForPsd = segment - mean(segment, 1, 'omitnan');
[psdMatrix, freqAxis] = pwelch(segmentForPsd, window, overlapSamples, nfft, Fs);

if plotInDb
    psdMatrix = 10 * log10(psdMatrix);
end

%% Plot combined PSD
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if showFigures
    figureVisibility = 'on';
else
    figureVisibility = 'off';
end

fig = figure('Color', 'w', 'Visible', figureVisibility, 'Position', [100 100 520 410]);
ax = axes('Parent', fig);
hold(ax, 'on');

colors = lines(size(psdMatrix, 2));
markers = {'o','d','s','^','v','>','<','p','h','x','+'};
markerEveryNSamples = max(1, floor(numel(freqAxis) / 25));

for ch = 1:size(psdMatrix, 2)
    marker = markers{mod(ch - 1, numel(markers)) + 1};
    plot(ax, freqAxis, psdMatrix(:, ch), ...
        'Color', colors(ch, :), ...
        'LineWidth', 1.4, ...
        'Marker', marker, ...
        'MarkerSize', 4.5, ...
        'MarkerIndices', 1:markerEveryNSamples:numel(freqAxis));
end

xlim(ax, frequencyRangeHz);
addFrequencyHighlight(ax, highlightFrequencyRangeHz, ...
    highlightBandColor, highlightBandAlpha, highlightLineColor);
grid(ax, 'off');
box(ax, 'off');
xlabel(ax, 'Frequency (Hz)', 'FontWeight', 'bold');

if plotInDb
    ylabel(ax, 'PSD (dB/Hz)', 'FontWeight', 'bold');
else
    ylabel(ax, 'PSD (\muV^2/Hz)', 'FontWeight', 'bold');
end

legend(ax, labels, 'Interpreter', 'none', 'Location', 'northeastoutside');
title(ax, sprintf('Combined PSD: %.2f to %.2f sec', actualStartTimeSec, actualEndTimeSec), ...
    'Interpreter', 'none');
set(ax, 'FontSize', 12, 'LineWidth', 1.2);

outputFilePath = buildOutputFilePath(filePath, 'Combined', 'PSD', outputFolder, outputImageFormat);
exportgraphics(fig, outputFilePath, 'Resolution', outputResolutionDpi);

%% Optional output in the base workspace
CombinedPSD_Result = struct();
CombinedPSD_Result.filePath = char(filePath);
CombinedPSD_Result.outputFilePath = outputFilePath;
CombinedPSD_Result.Fs = Fs;
CombinedPSD_Result.labels = labels;
CombinedPSD_Result.timePointSec = timePointSec;
CombinedPSD_Result.secondsBefore = secondsBefore;
CombinedPSD_Result.secondsAfter = secondsAfter;
CombinedPSD_Result.requestedStartTimeSec = requestedStartTimeSec;
CombinedPSD_Result.requestedEndTimeSec = requestedEndTimeSec;
CombinedPSD_Result.startTimeSec = actualStartTimeSec;
CombinedPSD_Result.endTimeSec = actualEndTimeSec;
CombinedPSD_Result.frequency = freqAxis;
CombinedPSD_Result.psd = psdMatrix;
assignin('base', 'CombinedPSD_Result', CombinedPSD_Result);

if closeFiguresAfterSave
    close(fig);
end

%% Local functions
function [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadEegFileForPsd(filePath, requestedLabels, requestedStartTimeSec, requestedEndTimeSec)
    [~, ~, ext] = fileparts(filePath);
    ext = lower(ext);

    switch ext
        case '.edf'
            [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadEdfForPsd( ...
                filePath, requestedLabels, requestedStartTimeSec, requestedEndTimeSec);
        case '.mat'
            [data, Fs, labels] = loadMatForPsd(filePath, requestedLabels);
            dataStartTimeSec = 0;
            recordDurationSec = size(data, 1) / Fs;
        otherwise
            error('Unsupported file extension "%s". Use an EDF or MAT file.', ext);
    end
end

function [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadEdfForPsd(filePath, requestedLabels, requestedStartTimeSec, requestedEndTimeSec)
    localEdfFolder = fullfile(fileparts(mfilename('fullpath')), 'PSG_chun_fft');
    if contains(path, localEdfFolder)
        rmpath(localEdfFolder);
    end

    if exist('edfinfo', 'file') ~= 2 || exist('edfread', 'file') ~= 2
        error('MATLAB edfinfo/edfread functions are required to read EDF files.');
    end

    info = edfinfo(filePath);
    allLabels = cellstr(info.SignalLabels);
    samplesPerRecord = double(info.NumSamples(:)');

    if isduration(info.DataRecordDuration)
        dataRecordDurationSec = seconds(info.DataRecordDuration);
    else
        dataRecordDurationSec = double(info.DataRecordDuration);
    end

    if isempty(requestedLabels)
        selectedIdx = find(~contains(lower(string(allLabels)), 'annotation'));
        if isempty(selectedIdx)
            selectedIdx = 1:numel(allLabels);
        end
    else
        selectedIdx = matchEdfLabels(allLabels, requestedLabels);
    end

    Fs = unique(samplesPerRecord(selectedIdx) ./ dataRecordDurationSec);
    if numel(Fs) ~= 1
        error(['Selected EDF channels have different sampling rates. ', ...
            'Choose channels with the same Fs or resample them before PSD.']);
    end

    Fs = double(Fs);
    labels = allLabels(selectedIdx);

    totalDataRecords = double(info.NumDataRecords);
    recordDurationSec = totalDataRecords * dataRecordDurationSec;
    firstDataRecord = max(1, floor(max(requestedStartTimeSec, 0) / dataRecordDurationSec) + 1);
    lastDataRecord = min(totalDataRecords, ceil(min(requestedEndTimeSec, recordDurationSec) / dataRecordDurationSec));

    if requestedEndTimeSec <= 0 || requestedStartTimeSec >= recordDurationSec || lastDataRecord < firstDataRecord
        error('Requested time period is outside the EDF recording duration %.3f sec.', recordDurationSec);
    end

    selectedDataRecords = firstDataRecord:lastDataRecord;
    dataStartTimeSec = (firstDataRecord - 1) * dataRecordDurationSec;

    edfTable = edfread(filePath, ...
        'SelectedSignals', string(labels), ...
        'SelectedDataRecords', selectedDataRecords);
    data = timetableToMatrixForPsd(edfTable);
end

function data = timetableToMatrixForPsd(edfTable)
    nChannels = width(edfTable);
    channelVectors = cell(1, nChannels);
    nSamples = [];

    for ch = 1:nChannels
        variableName = edfTable.Properties.VariableNames{ch};
        recordCells = edfTable.(variableName);

        if iscell(recordCells)
            channelData = vertcat(recordCells{:});
        else
            channelData = recordCells;
        end

        channelData = double(channelData(:));

        if isempty(nSamples)
            nSamples = numel(channelData);
        elseif numel(channelData) ~= nSamples
            error(['Selected EDF channels do not contain the same number of samples. ', ...
                'Choose channels with the same sampling rate and duration.']);
        end

        channelVectors{ch} = channelData;
    end

    data = zeros(nSamples, nChannels);
    for ch = 1:nChannels
        data(:, ch) = channelVectors{ch};
    end
end

function idx = matchEdfLabels(labels, requestedLabels)
    if isnumeric(requestedLabels)
        idx = requestedLabels(:)';
        if any(idx < 1) || any(idx > numel(labels)) || any(fix(idx) ~= idx)
            error('EDF channel indices must be integers from 1 to %d.', numel(labels));
        end
        return;
    end

    labelsString = string(labels);
    cleanLabels = regexprep(labelsString, '\W', '');
    requestedString = string(requestedLabels);
    cleanRequested = regexprep(requestedString, '\W', '');
    idx = zeros(1, numel(requestedString));

    for i = 1:numel(requestedString)
        match = find(labelsString == requestedString(i), 1);

        if isempty(match)
            match = find(lower(labelsString) == lower(requestedString(i)), 1);
        end

        if isempty(match)
            match = find(cleanLabels == cleanRequested(i), 1);
        end

        if isempty(match)
            match = find(lower(cleanLabels) == lower(cleanRequested(i)), 1);
        end

        if isempty(match)
            error('EDF channel "%s" was not found.', requestedString(i));
        end

        idx(i) = match;
    end
end

function [data, Fs, labels] = loadMatForPsd(filePath, requestedLabels)
    S = load(filePath);
    Fs = findSamplingRate(S);

    if isfield(S, 'EEG') && isstruct(S.EEG) && isfield(S.EEG, 'data')
        [data, labels] = unpackEeglabStruct(S.EEG, requestedLabels);
        return;
    end

    if ~isempty(requestedLabels)
        [data, labels] = loadRequestedMatChannels(S, requestedLabels);
        return;
    end

    [data, labels] = loadDetectedMatChannels(S);
end

function Fs = findSamplingRate(S)
    fsNames = {'Fs', 'fs', 'srate', 'sampleRate', 'samplingRate', 'SamplingRate'};

    for i = 1:numel(fsNames)
        name = fsNames{i};
        if isfield(S, name) && isnumeric(S.(name)) && isscalar(S.(name))
            Fs = double(S.(name));
            return;
        end
    end

    if isfield(S, 'EEG') && isstruct(S.EEG) && isfield(S.EEG, 'srate')
        Fs = double(S.EEG.srate);
        return;
    end

    error(['Could not find sampling rate in MAT file. Add a scalar variable ', ...
        'named Fs, fs, srate, sampleRate, or samplingRate.']);
end

function [data, labels] = unpackEeglabStruct(EEG, requestedLabels)
    rawData = double(EEG.data);

    if ~ismatrix(rawData)
        rawData = reshape(rawData, size(rawData, 1), []);
    end

    if isfield(EEG, 'chanlocs') && isfield(EEG.chanlocs, 'labels')
        labels = {EEG.chanlocs.labels};
    else
        labels = arrayfun(@(x) sprintf('Ch%d', x), 1:size(rawData, 1), ...
            'UniformOutput', false);
    end

    data = rawData';

    if ~isempty(requestedLabels)
        idx = matchLabels(labels, requestedLabels);
        data = data(:, idx);
        labels = labels(idx);
    end
end

function [data, labels] = loadRequestedMatChannels(S, requestedLabels)
    labels = cellstr(requestedLabels);
    data = [];

    for i = 1:numel(labels)
        name = labels{i};
        matchName = findMatchingField(S, name);

        if strlength(matchName) == 0
            error('Channel "%s" was not found as a MAT variable.', name);
        end

        matchName = char(matchName);
        x = S.(matchName);
        if ~isnumeric(x) || ~isvector(x)
            error('MAT variable "%s" must be a numeric vector channel.', matchName);
        end

        data(:, i) = double(x(:)); %#ok<AGROW>
        labels{i} = char(matchName);
    end
end

function [data, labels] = loadDetectedMatChannels(S)
    names = fieldnames(S);
    vectorNames = {};
    vectorLengths = [];

    for i = 1:numel(names)
        name = names{i};
        value = S.(name);

        if isnumeric(value) && isvector(value) && numel(value) > 1
            vectorNames{end + 1} = name; %#ok<AGROW>
            vectorLengths(end + 1) = numel(value); %#ok<AGROW>
        end
    end

    if ~isempty(vectorNames)
        targetLength = max(vectorLengths);
        keep = vectorLengths == targetLength;
        labels = vectorNames(keep);
        data = zeros(targetLength, numel(labels));

        for i = 1:numel(labels)
            data(:, i) = double(S.(labels{i})(:));
        end

        return;
    end

    if isfield(S, 'data') && isnumeric(S.data) && ismatrix(S.data)
        data = double(S.data);
        if size(data, 1) < size(data, 2)
            data = data';
        end

        labels = arrayfun(@(x) sprintf('Ch%d', x), 1:size(data, 2), ...
            'UniformOutput', false);
        return;
    end

    error(['No channel data was detected in the MAT file. Expected named ', ...
        'numeric channel vectors, an EEG.data struct, or a numeric data matrix.']);
end

function fieldName = findMatchingField(S, requestedName)
    names = string(fieldnames(S));
    requestedName = string(requestedName);

    exact = find(names == requestedName, 1);
    if ~isempty(exact)
        fieldName = names(exact);
        return;
    end

    caseInsensitive = find(lower(names) == lower(requestedName), 1);
    if ~isempty(caseInsensitive)
        fieldName = names(caseInsensitive);
        return;
    end

    fieldName = "";
end

function idx = matchLabels(labels, requestedLabels)
    labelsString = string(labels);
    requestedString = string(requestedLabels);
    idx = zeros(1, numel(requestedString));

    for i = 1:numel(requestedString)
        exact = find(labelsString == requestedString(i), 1);

        if isempty(exact)
            exact = find(lower(labelsString) == lower(requestedString(i)), 1);
        end

        if isempty(exact)
            error('Channel "%s" was not found.', requestedString(i));
        end

        idx(i) = exact;
    end
end

function addFrequencyHighlight(ax, frequencyRangeHz, bandColor, bandAlpha, lineColor)
    if isempty(frequencyRangeHz)
        return;
    end

    if numel(frequencyRangeHz) ~= 2 || frequencyRangeHz(2) <= frequencyRangeHz(1)
        error('highlightFrequencyRangeHz must be [] or a two-element range like [25 50].');
    end

    xLimits = xlim(ax);
    bandStart = max(frequencyRangeHz(1), xLimits(1));
    bandEnd = min(frequencyRangeHz(2), xLimits(2));

    if bandEnd <= bandStart
        return;
    end

    yLimits = ylim(ax);
    wasHolding = ishold(ax);
    hold(ax, 'on');

    bandPatch = patch(ax, ...
        [bandStart bandEnd bandEnd bandStart], ...
        [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
        bandColor, ...
        'FaceAlpha', bandAlpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    uistack(bandPatch, 'bottom');

    xline(ax, bandStart, '--', 'Color', lineColor, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
    xline(ax, bandEnd, '--', 'Color', lineColor, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
    ylim(ax, yLimits);

    if ~wasHolding
        hold(ax, 'off');
    end
end

function outputFilePath = buildOutputFilePath(filePath, channelLabel, suffix, outputFolder, outputImageFormat)
    [~, eegFileName] = fileparts(filePath);
    outputFileName = sprintf('%s_%s_%s.%s', ...
        sanitizeFileName(eegFileName), ...
        sanitizeFileName(channelLabel), ...
        suffix, ...
        sanitizeFileName(outputImageFormat));
    outputFilePath = fullfile(outputFolder, outputFileName);
end

function safeName = sanitizeFileName(value)
    safeName = char(string(value));
    safeName = regexprep(safeName, '[<>:"/\\|?*]', '_');
    safeName = regexprep(safeName, '\s+', '_');
    safeName = regexprep(safeName, '_+', '_');
    safeName = strtrim(safeName);

    if isempty(safeName)
        safeName = 'unnamed';
    end
end
