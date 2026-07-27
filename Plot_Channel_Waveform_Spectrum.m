%% Plot original waveform and power spectrum for selected EEG channels
% Edit the settings in this section, then run the script.
% EDF files are read with MATLAB edfinfo/edfread.

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
    filePath = "";                 % Example: "D:\EEG Files\Yeh\0413_1551-1601.EDF"
end
if ~exist('channelNames', 'var')
    if exist('channelName', 'var')
        channelNames = channelName;        % Backward-compatible with the old setting name
    else
        channelNames = ["F3" "F7"];       % Example: ["F3", "F4", "O1"]
        % channelNames = ["F3-AVG" "F4-AVG"];
        % channelNames = ["Fpz"];
    end
end
if ~exist('channelIndex', 'var')
    channelIndex = [];              % Example: [3 4]; used only when channelNames is ""
end
if ~exist('startTimeSec', 'var')
    startTimeSec = 26441.5;               % Segment start time in seconds
end
if ~exist('endTimeSec', 'var')
    endTimeSec = 26446.5;                % Segment end time in seconds
end
if ~exist('timePointSec', 'var')
    timePointSec = [];                 % Example: 100; [] uses startTimeSec/endTimeSec instead
end
if ~exist('secondsBefore', 'var')
    secondsBefore = 2.5;                % Seconds plotted before timePointSec
end
if ~exist('secondsAfter', 'var')
    secondsAfter = 2.5;                 % Seconds plotted after timePointSec
end

if ~exist('frequencyRangeHz', 'var')
    frequencyRangeHz = [20 50];     % Frequency range shown in the spectrum plot
end
if ~exist('spectrogramWindowSec', 'var')
    spectrogramWindowSec = 1;     % Shorter windows improve time resolution
end
if ~exist('spectrogramOverlapFraction', 'var')
    spectrogramOverlapFraction = 0.5;
end
if ~exist('nfft', 'var')
    nfft = [];                      % [] chooses an automatic value
end
if ~exist('minMaxNormalizeSpectrum', 'var')
    minMaxNormalizeSpectrum = true; % true gives a 0 to 1 color scale like the example
end
if ~exist('showPowerDb', 'var')
    showPowerDb = false;            % used only when minMaxNormalizeSpectrum is false
end
if ~exist('timeTickFormat', 'var')
    timeTickFormat = '%.3f';         % Plain x-axis time labels, no scientific notation
end
if ~exist('outputFolder', 'var')
    outputFolder = fullfile(pwd, 'Yeh_PSDandSpectrogram/2大/');
    if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
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
% if ~exist('waveformYLim', 'var')
%     waveformYLim = [-500 500];
% end
%% Select file
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

channelNames = string(channelNames);
usingChannelNames = any(strlength(channelNames) > 0);

if ~usingChannelNames
    if isempty(channelIndex)
        channelIndex = 1;
    end
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if showFigures
    figureVisibility = 'on';
else
    figureVisibility = 'off';
end

WaveformSpectrum_Result = struct();

%% Determine requested time period
if isempty(timePointSec)
    requestedStartTimeSec = startTimeSec;
    requestedEndTimeSec = endTimeSec;
else
    if ~isscalar(timePointSec)
        error('timePointSec must be a scalar value in seconds.');
    end
    timePointSec = double(timePointSec);
    if timePointSec < 0 || timePointSec > recordDurationSec
        error('timePointSec %.3f is outside the recording duration %.3f sec.', ...
            timePointSec, recordDurationSec);
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
[eegData, Fs, labels, dataStartTimeSec, recordDurationSec] = loadSelectedEegChannels( ...
    filePath, channelNames, channelIndex, usingChannelNames, ...
    requestedStartTimeSec, requestedEndTimeSec);
nChannelsToPlot = numel(labels);

%% Select time period once
nSamples = size(eegData, 1);

if requestedStartTimeSec < 0
    warning('Requested start time %.3f sec is before the recording start. Plotting from 0 sec.', ...
        requestedStartTimeSec);
end

if requestedEndTimeSec > recordDurationSec
    warning(['Requested end time %.3f sec is after the recording end %.3f sec. ', ...
        'Plotting until the end of the recording.'], ...
        requestedEndTimeSec, recordDurationSec);
end

if requestedEndTimeSec <= 0 || requestedStartTimeSec >= recordDurationSec || ...
        requestedEndTimeSec <= requestedStartTimeSec
    error('Use a valid time period: 0 <= startTimeSec < endTimeSec.');
end

if spectrogramOverlapFraction < 0 || spectrogramOverlapFraction >= 1
    error('spectrogramOverlapFraction must be >= 0 and < 1.');
end
if numel(frequencyRangeHz) ~= 2 || frequencyRangeHz(1) < 0 || frequencyRangeHz(2) <= frequencyRangeHz(1)
    error('frequencyRangeHz must be a two-element range like [20 50].');
end

firstSample = max(1, floor((requestedStartTimeSec - dataStartTimeSec) * Fs) + 1);
lastSample = min(nSamples, floor((requestedEndTimeSec - dataStartTimeSec) * Fs));
segments = eegData(firstSample:lastSample, :);
actualStartTimeSec = dataStartTimeSec + (firstSample - 1) / Fs;
actualEndTimeSec = dataStartTimeSec + lastSample / Fs;
timeAxis = actualStartTimeSec + (0:size(segments, 1) - 1) ./ Fs;

if size(segments, 1) < 4
    error('Selected segment is too short for spectrum calculation.');
end

%% Shared spectrogram parameters
windowLengthSamples = max(4, round(spectrogramWindowSec * Fs));
windowLengthSamples = min(windowLengthSamples, size(segments, 1));
overlapSamples = floor(windowLengthSamples * spectrogramOverlapFraction);

if isempty(nfft)
    effectiveNfft = max(256, 2 ^ nextpow2(windowLengthSamples));
else
    effectiveNfft = nfft;
end

window = hamming(windowLengthSamples, 'periodic');

%% Compute and plot each channel
for ch = 1:nChannelsToPlot
    label = labels{ch};
    segment = segments(:, ch);
    
    % bandpassRangeHz = [25 50];   % edit as needed
    % segmentForSpectrum = bandpass(segment, bandpassRangeHz, Fs);
    % segmentForSpectrum = segmentForSpectrum - mean(segmentForSpectrum, 'omitnan');
    segmentForSpectrum = segment - mean(segment, 'omitnan');

    [~, frequencyAxis, spectrumTimeAxis, psdSpectrum] = spectrogram( ...
        segmentForSpectrum, ...
        window, ...
        overlapSamples, ...
        effectiveNfft, ...
        Fs, ...
        'psd');
    
    spectrumTimeAxis = spectrumTimeAxis + actualStartTimeSec;
    
    freqMask = frequencyAxis >= frequencyRangeHz(1) & frequencyAxis <= frequencyRangeHz(2);
    plotFrequencyAxis = frequencyAxis(freqMask);
    rawPlotPowerSpectrum = psdSpectrum(freqMask, :);
    [plotPowerSpectrum, powerSpectrum, spectrumClim, climLabel] = formatSpectrumForDisplay( ...
        rawPlotPowerSpectrum, minMaxNormalizeSpectrum, showPowerDb);

    fig = figure('Color', 'w', 'Visible', figureVisibility, 'Position', [100 100 980 420]);
    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    axWave = nexttile;
    plot(axWave, timeAxis, segment, 'k', 'LineWidth', 0.8);
    title(axWave, label, 'Interpreter', 'none', 'FontSize', 18, 'FontWeight', 'bold');
    xlim(axWave, [timeAxis(1), timeAxis(end)]);
    formatPlainTimeAxis(axWave, timeTickFormat);
    % if ~isempty(waveformYLim)
    %     ylim(axWave, waveformYLim);
    % end
    ylabel(axWave, 'Amplitude (\muV)');
    grid(axWave, 'on');
    box(axWave, 'on');
    if ~isempty(timePointSec)
        xline(axWave, timePointSec, 'r--', 'LineWidth', 1.2);
    end
    
    axSpectrum = nexttile;
    imagesc(axSpectrum, spectrumTimeAxis, plotFrequencyAxis, plotPowerSpectrum);
    axis(axSpectrum, 'xy');
    colormap(axSpectrum, jet);
    colorbar(axSpectrum);
    clim(axSpectrum, spectrumClim);
    xlim(axSpectrum, [timeAxis(1), timeAxis(end)]);
    formatPlainTimeAxis(axSpectrum, timeTickFormat);
    ylim(axSpectrum, frequencyRangeHz);
    xlabel(axSpectrum, 'Time (s)');
    ylabel(axSpectrum, 'Frequency (Hz)');
    box(axSpectrum, 'on');
    if ~isempty(timePointSec)
        xline(axSpectrum, timePointSec, 'r--', 'LineWidth', 1.2);
    end

    title(axSpectrum, climLabel, 'FontWeight', 'normal');

    outputFilePath = buildOutputFilePath(filePath, label, 'Waveform_Spectrogram', ...
        outputFolder, outputImageFormat);
    exportgraphics(fig, outputFilePath, 'Resolution', outputResolutionDpi);

    WaveformSpectrum_Result(ch).filePath = char(filePath);
    WaveformSpectrum_Result(ch).outputFilePath = outputFilePath;
    WaveformSpectrum_Result(ch).channel = label;
    WaveformSpectrum_Result(ch).Fs = Fs;
    WaveformSpectrum_Result(ch).timePointSec = timePointSec;
    WaveformSpectrum_Result(ch).secondsBefore = secondsBefore;
    WaveformSpectrum_Result(ch).secondsAfter = secondsAfter;
    WaveformSpectrum_Result(ch).requestedStartTimeSec = requestedStartTimeSec;
    WaveformSpectrum_Result(ch).requestedEndTimeSec = requestedEndTimeSec;
    WaveformSpectrum_Result(ch).startTimeSec = actualStartTimeSec;
    WaveformSpectrum_Result(ch).endTimeSec = actualEndTimeSec;
    WaveformSpectrum_Result(ch).time = timeAxis(:);
    WaveformSpectrum_Result(ch).waveform = segment(:);
    WaveformSpectrum_Result(ch).spectrumTime = spectrumTimeAxis(:);
    WaveformSpectrum_Result(ch).frequency = plotFrequencyAxis(:);
    WaveformSpectrum_Result(ch).powerSpectrum = powerSpectrum;

    if closeFiguresAfterSave
        close(fig);
    end
end

%% Optional output in the base workspace
assignin('base', 'WaveformSpectrum_Result', WaveformSpectrum_Result);

%% Local functions
function [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadSelectedEegChannels(filePath, channelNames, channelIndex, usingChannelNames, requestedStartTimeSec, requestedEndTimeSec)
    [~, ~, ext] = fileparts(filePath);
    ext = lower(ext);

    switch ext
        case '.edf'
            [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadSelectedEdfChannels( ...
                filePath, channelNames, channelIndex, usingChannelNames, ...
                requestedStartTimeSec, requestedEndTimeSec);
        case '.mat'
            [data, Fs, labels] = loadSelectedMatChannels(filePath, channelNames, channelIndex, usingChannelNames);
            dataStartTimeSec = 0;
            recordDurationSec = size(data, 1) / Fs;
        otherwise
            error('Unsupported file extension "%s". Use an EDF or MAT file.', ext);
    end
end

function [data, Fs, labels, dataStartTimeSec, recordDurationSec] = loadSelectedEdfChannels(filePath, channelNames, channelIndex, usingChannelNames, requestedStartTimeSec, requestedEndTimeSec)
    localEdfFolder = fullfile(fileparts(mfilename('fullpath')), 'PSG_chun_fft');
    if contains(path, localEdfFolder)
        rmpath(localEdfFolder);
    end

    if exist('edfinfo', 'file') ~= 2 || exist('edfread', 'file') ~= 2
        error('MATLAB edfinfo/edfread functions are required to read EDF files.');
    end

    info = edfinfo(filePath);
    allLabels = cellstr(info.SignalLabels);

    if usingChannelNames
        selectedIdx = resolveChannelIndices(allLabels, channelNames);
    else
        selectedIdx = channelIndex(:)';
    end

    if any(selectedIdx < 1) || any(selectedIdx > numel(allLabels)) || any(fix(selectedIdx) ~= selectedIdx)
        error('channelIndex values must be integers from 1 to %d.', numel(allLabels));
    end

    samplesPerRecord = double(info.NumSamples(selectedIdx));
    if isduration(info.DataRecordDuration)
        dataRecordDurationSec = seconds(info.DataRecordDuration);
    else
        dataRecordDurationSec = double(info.DataRecordDuration);
    end

    FsValues = unique(samplesPerRecord ./ dataRecordDurationSec);
    if numel(FsValues) ~= 1
        error(['Selected EDF channels have different sampling rates. ', ...
            'Choose channels with the same Fs or resample them before plotting.']);
    end

    Fs = double(FsValues);
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
    data = timetableToMatrixForSpectrum(edfTable);
end

function data = timetableToMatrixForSpectrum(edfTable)
    nChannels = width(edfTable);
    channelVectors = cell(1, nChannels);
    nSamples = [];

    for ch = 1:nChannels
        channelData = timetableColumnToVector(edfTable, ch);

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

function [data, Fs, labels] = loadSelectedMatChannels(filePath, channelNames, channelIndex, usingChannelNames)
    S = load(filePath);
    Fs = findSamplingRate(S);

    if isfield(S, 'EEG') && isstruct(S.EEG) && isfield(S.EEG, 'data')
        [allData, allLabels] = unpackEeglabStruct(S.EEG);
    else
        [allData, allLabels] = unpackMatChannels(S);
    end

    if usingChannelNames
        selectedIdx = resolveChannelIndices(allLabels, channelNames);
    else
        selectedIdx = channelIndex(:)';
    end

    if any(selectedIdx < 1) || any(selectedIdx > numel(allLabels)) || any(fix(selectedIdx) ~= selectedIdx)
        error('channelIndex values must be integers from 1 to %d.', numel(allLabels));
    end

    data = allData(:, selectedIdx);
    labels = allLabels(selectedIdx);
end

function idx = resolveChannelIndices(labels, requestedLabels)
    requestedLabels = string(requestedLabels);
    idx = zeros(1, numel(requestedLabels));

    for i = 1:numel(requestedLabels)
        idx(i) = matchOneLabel(labels, requestedLabels(i));
    end
end

function signal = timetableColumnToVector(edfTable, columnIndex)
    variableName = edfTable.Properties.VariableNames{columnIndex};
    recordCells = edfTable.(variableName);

    if iscell(recordCells)
        signal = vertcat(recordCells{:});
    else
        signal = recordCells;
    end

    signal = double(signal(:));
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

function [data, labels] = unpackEeglabStruct(EEG)
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
end

function [data, labels] = unpackMatChannels(S)
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

function idx = matchOneLabel(labels, requestedLabel)
    labelsString = string(labels);
    requestedString = string(requestedLabel);
    cleanLabels = regexprep(labelsString, '\W', '');
    cleanRequested = regexprep(requestedString, '\W', '');

    idx = find(labelsString == requestedString, 1);

    if isempty(idx)
        idx = find(lower(labelsString) == lower(requestedString), 1);
    end

    if isempty(idx)
        idx = find(cleanLabels == cleanRequested, 1);
    end

    if isempty(idx)
        idx = find(lower(cleanLabels) == lower(cleanRequested), 1);
    end

    if isempty(idx)
        error('Channel "%s" was not found.', requestedString);
    end
end

function [plotPowerSpectrum, storedPowerSpectrum, spectrumClim, climLabel] = formatSpectrumForDisplay(rawPowerSpectrum, minMaxNormalizeSpectrum, showPowerDb)
    plotPowerSpectrum = rawPowerSpectrum;

    if minMaxNormalizeSpectrum
        minPower = min(plotPowerSpectrum(:));
        maxPower = max(plotPowerSpectrum(:));

        if maxPower > minPower
            plotPowerSpectrum = (plotPowerSpectrum - minPower) ./ (maxPower - minPower);
        else
            plotPowerSpectrum(:) = 0;
        end

        storedPowerSpectrum = plotPowerSpectrum;
        spectrumClim = [0 1];
        climLabel = 'Normalized power';
        return;
    end

    if showPowerDb
        plotPowerSpectrum = 10 * log10(plotPowerSpectrum + eps);
        climLabel = 'Power (dB/Hz)';
    else
        climLabel = 'PSD (\muV^2/Hz)';
    end

    storedPowerSpectrum = plotPowerSpectrum;
    spectrumClim = percentileColorLimits(plotPowerSpectrum, 5, 99);
end

function limits = percentileColorLimits(values, lowPercentile, highPercentile)
    finiteValues = values(isfinite(values));

    if isempty(finiteValues)
        limits = [0 1];
        return;
    end

    limits = prctile(finiteValues(:), [lowPercentile highPercentile]);

    if limits(1) == limits(2)
        limits = limits + [-1 1] * max(eps, abs(limits(1)) * 0.01);
    end
end

function formatPlainTimeAxis(ax, tickFormat)
    ax.XAxis.Exponent = 0;
    ax.XTickLabel = compose(tickFormat, ax.XTick);
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
