%% Plot whole-span and time-resolved gamma-band scalp topographies
% First run Plot_Combined_Channel_PSD.m to load the selected channels and
% five-second segment. This unified script can create:
%   1. one topoplot for the entire selected span,
%   2. one topoplot per non-overlapping epoch, or
%   3. both figures.

clc;
clearvars -except CombinedPSD_Result commonMapLimits plotMode gammaBandHz analysisDurationSec ...
    epochLengthSec displayMode baselineEpochIndices channelLocationFile ...
    outputFolder outputResolutionDpi showFigures closeFiguresAfterSave EEG;

%% User settings
if ~exist('gammaBandHz', 'var')
    gammaBandHz = [25 50];
end
if ~exist('analysisDurationSec', 'var')
    analysisDurationSec = 5;
end
if ~exist('epochLengthSec', 'var')
    epochLengthSec = 0.5;
end
if ~exist('plotMode', 'var')
    % "both", "whole-span", or "epochs"
    plotMode = "both";
end
if ~exist('displayMode', 'var')
    % "absolute-linear" shows raw integrated gamma power.
    % "absolute-db" shows gamma power in dB.
    % "relative-db" shows dB change relative to baselineEpochIndices.
    displayMode = "absolute-db";
end
if ~exist('baselineEpochIndices', 'var')
    baselineEpochIndices = 1:2; % First second; used only for "relative-db"
end
if ~exist('channelLocationFile', 'var')
    % Optional EEGLAB-readable .ced, .locs, .sfp, or similar file.
    % Leave empty to use EEG.chanlocs from a source MAT file, when
    % available, or EEGLAB's standard channel-location lookup.
    channelLocationFile = "";
end
if ~exist('outputFolder', 'var')
    outputFolder = "";
end
if ~exist('outputResolutionDpi', 'var')
    outputResolutionDpi = 300;
end
if ~exist('showFigures', 'var')
    showFigures = true;
end
if ~exist('closeFiguresAfterSave', 'var')
    closeFiguresAfterSave = false;
end

%% Validate input
if ~exist('CombinedPSD_Result', 'var') || ~isstruct(CombinedPSD_Result)
    error(['CombinedPSD_Result was not found. Run ', ...
        'Plot_Combined_Channel_PSD.m first.']);
end

requiredFields = {'segment', 'Fs', 'labels', 'filePath', 'startTimeSec'};
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};
    if ~isfield(CombinedPSD_Result, fieldName)
        error(['CombinedPSD_Result.%s was not found. Run the updated ', ...
            'Plot_Combined_Channel_PSD.m first.'], fieldName);
    end
end

segment = double(CombinedPSD_Result.segment);
Fs = double(CombinedPSD_Result.Fs);
labels = cellstr(string(CombinedPSD_Result.labels));
segmentStartTimeSec = double(CombinedPSD_Result.startTimeSec);

if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('CombinedPSD_Result.Fs must be a positive sampling rate.');
end
if size(segment, 2) ~= numel(labels)
    error('The segment channel count does not match the number of labels.');
end
if numel(gammaBandHz) ~= 2 || gammaBandHz(2) <= gammaBandHz(1)
    error('gammaBandHz must be an increasing range such as [30 50].');
end
if gammaBandHz(1) < 0 || gammaBandHz(2) > Fs / 2
    error('gammaBandHz must lie between 0 and the Nyquist frequency %.2f Hz.', Fs / 2);
end
if analysisDurationSec <= 0 || epochLengthSec <= 0
    error('analysisDurationSec and epochLengthSec must be positive.');
end

analysisSamples = round(analysisDurationSec * Fs);
epochSamples = round(epochLengthSec * Fs);
if analysisSamples > size(segment, 1)
    error(['The saved segment contains %.3f seconds, but %.3f seconds were ', ...
        'requested. Increase the period in Plot_Combined_Channel_PSD.m.'], ...
        size(segment, 1) / Fs, analysisDurationSec);
end
if epochSamples < 4
    error('Each epoch must contain at least four samples.');
end

numberOfEpochsExact = analysisSamples / epochSamples;
numberOfEpochs = round(numberOfEpochsExact);
if abs(numberOfEpochsExact - numberOfEpochs) > 1e-10
    error(['analysisDurationSec must contain a whole number of epochs. ', ...
        'For example, 5 seconds / 0.5 seconds = 10 epochs.']);
end

analysisSegment = segment(1:analysisSamples, :);

%% Calculate gamma power separately for every epoch
if exist('bandpower', 'file') ~= 2 || exist('pwelch', 'file') ~= 2
    error('MATLAB Signal Processing Toolbox is required.');
end

nfft = max(256, 2 ^ nextpow2(epochSamples));
epochWindow = hamming(epochSamples, 'periodic');
gammaPower = zeros(numberOfEpochs, size(analysisSegment, 2));
epochStartTimeSec = zeros(numberOfEpochs, 1);
epochEndTimeSec = zeros(numberOfEpochs, 1);
epochCenterTimeSec = zeros(numberOfEpochs, 1);

for epochIndex = 1:numberOfEpochs
    firstSample = (epochIndex - 1) * epochSamples + 1;
    lastSample = firstSample + epochSamples - 1;
    epochData = analysisSegment(firstSample:lastSample, :);
    epochData = detrend(epochData, 0);

    % A single Hamming-windowed PSD is used because the entire epoch is
    % only 0.5 seconds. Subdividing it further would reduce true frequency
    % resolution below what is useful for the 30-50 Hz band.
    [epochPsd, frequency] = pwelch( ...
        epochData, epochWindow, 0, nfft, Fs);
    gammaPower(epochIndex, :) = bandpower( ...
        epochPsd, frequency, gammaBandHz, "psd");

    epochStartTimeSec(epochIndex) = ...
        segmentStartTimeSec + (firstSample - 1) / Fs;
    epochEndTimeSec(epochIndex) = ...
        segmentStartTimeSec + lastSample / Fs;
    epochCenterTimeSec(epochIndex) = ...
        (epochStartTimeSec(epochIndex) + epochEndTimeSec(epochIndex)) / 2;
end

%% Select absolute or baseline-relative display values
displayMode = lower(string(displayMode));
switch displayMode
    case "absolute-linear"
        topoplotValues = gammaPower;
        colorbarLabel = 'Gamma power (\muV^2)';

    case "absolute-db"
        topoplotValues = 10 * log10(max(gammaPower, realmin('double')));
        colorbarLabel = 'Gamma power (dB re 1 \muV^2)';

    case "relative-db"
        if isempty(baselineEpochIndices) || ...
                any(baselineEpochIndices < 1) || ...
                any(baselineEpochIndices > numberOfEpochs) || ...
                any(fix(baselineEpochIndices) ~= baselineEpochIndices)
            error('baselineEpochIndices must contain valid integer epoch indices.');
        end

        baselinePower = mean(gammaPower(baselineEpochIndices, :), 1);
        topoplotValues = 10 * log10( ...
            max(gammaPower, realmin('double')) ./ ...
            max(baselinePower, realmin('double')));
        colorbarLabel = 'Gamma change from baseline (dB)';

    otherwise
        error(['displayMode must be "absolute-linear", "absolute-db", ', ...
            'or "relative-db".']);
end

%% Match power values with EEGLAB channel locations
if exist('topoplot', 'file') ~= 2 || exist('pop_chanedit', 'file') ~= 2
    error(['EEGLAB was not found on the MATLAB path. Start EEGLAB or add ', ...
        'its folder to the MATLAB path, then run this script again.']);
end

chanlocs = resolveTimeTopoplotLocations( ...
    labels, CombinedPSD_Result.filePath, channelLocationFile);

%% Select which figures to create
plotMode = lower(string(plotMode));
validPlotModes = ["both", "whole-span", "epochs"];
if ~isscalar(plotMode) || ~any(plotMode == validPlotModes)
    error('plotMode must be "both", "whole-span", or "epochs".');
end

plotWholeSpan = plotMode == "both" || plotMode == "whole-span";
plotEpochs = plotMode == "both" || plotMode == "epochs";

if showFigures
    figureVisibility = 'on';
else
    figureVisibility = 'off';
end

%% Prepare output paths
if strlength(string(outputFolder)) == 0
    if isfield(CombinedPSD_Result, 'outputFilePath')
        outputFolder = fileparts(CombinedPSD_Result.outputFilePath);
    else
        outputFolder = fullfile(pwd, 'Gamma_Topoplot');
    end
end
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

[~, sourceName] = fileparts(CombinedPSD_Result.filePath);
safeSourceName = sanitizeTimeTopoplotFileName(sourceName);
wholeSpanOutputFileName = sprintf( ...
    '%s_Gamma_%g-%gHz_WholeSpan_Topoplot.png', ...
    safeSourceName, gammaBandHz(1), gammaBandHz(2));
epochOutputFileName = sprintf( ...
    '%s_Gamma_%g-%gHz_%gsEpoch_Topoplots.png', ...
    safeSourceName, ...
    gammaBandHz(1), ...
    gammaBandHz(2), ...
    epochLengthSec);
wholeSpanOutputFilePath = fullfile(outputFolder, wholeSpanOutputFileName);
epochOutputFilePath = fullfile(outputFolder, epochOutputFileName);

%% Plot one map for the whole selected span
wholeSpanGammaPower = mean(gammaPower, 1);
if displayMode == "relative-db"
    wholeSpanTopoplotValues = 10 * log10( ...
        max(wholeSpanGammaPower, realmin('double')) ./ ...
        max(baselinePower, realmin('double')));
elseif displayMode == "absolute-linear"
    wholeSpanTopoplotValues = wholeSpanGammaPower;
else
    wholeSpanTopoplotValues = 10 * log10( ...
        max(wholeSpanGammaPower, realmin('double')));
end

% Use one color scale for every map in this section. Include both the
% whole-span and epoch values so that equal colors represent equal power
% without clipping any epoch/channel values.
allSectionTopoplotValues = [ ...
    topoplotValues(:); ...
    wholeSpanTopoplotValues(:)];
if displayMode == "relative-db"
    maximumAbsoluteValue = max(abs(allSectionTopoplotValues));
    if ~exist('commonMapLimits', 'var')
        commonMapLimits = [-maximumAbsoluteValue maximumAbsoluteValue];
    end
else
    if ~exist('commonMapLimits', 'var')
        commonMapLimits = [ ...
            min(allSectionTopoplotValues) max(allSectionTopoplotValues)];
    end
end
if commonMapLimits(1) == commonMapLimits(2)
    limitPadding = max(abs(commonMapLimits(1)) * 0.01, 1e-12);
    commonMapLimits = commonMapLimits + [-limitPadding limitPadding];
end

mapLimits = commonMapLimits;
wholeSpanMapLimits = commonMapLimits;
numberOfSaturatedEpochValues = 0;

wholeSpanFigure = [];
if plotWholeSpan
    wholeSpanFigure = figure( ...
        'Color', 'w', ...
        'Visible', figureVisibility, ...
        'Position', [100 100 600 520]);
    wholeSpanAxes = axes('Parent', wholeSpanFigure);
    axes(wholeSpanAxes); 
    topoplot(wholeSpanTopoplotValues(:), chanlocs, ...
        'electrodes', 'labels', ...
        'style', 'both');
    clim(wholeSpanAxes, wholeSpanMapLimits);

    colormap(wholeSpanFigure, jet(256));

    wholeSpanColorbar = colorbar(wholeSpanAxes);
    wholeSpanColorbar.Label.String = colorbarLabel;
    title(wholeSpanAxes, sprintf( ...
        ['Whole-span Gamma Topography (%.1f-%.1f Hz)\n', ...
        '%.2f-%.2f s'], ...
        gammaBandHz(1), gammaBandHz(2), ...
        epochStartTimeSec(1), epochEndTimeSec(end)));
    exportgraphics( ...
        wholeSpanFigure, wholeSpanOutputFilePath, ...
        'Resolution', outputResolutionDpi);
end

%% Plot every epoch with one common color scale
epochFigure = [];
if plotEpochs
    numberOfColumns = min(5, numberOfEpochs);
    numberOfRows = ceil(numberOfEpochs / numberOfColumns);
    figureWidth = max(900, 190 * numberOfColumns);
    figureHeight = max(480, 230 * numberOfRows);

    epochFigure = figure( ...
        'Color', 'w', ...
        'Visible', figureVisibility, ...
        'Position', [50 50 figureWidth figureHeight]);

    lastAxes = [];
    for epochIndex = 1:numberOfEpochs
        lastAxes = subplot( ...
            numberOfRows, numberOfColumns, epochIndex, ...
            'Parent', epochFigure);
        axes(lastAxes); %#ok<LAXES>
        topoplot(topoplotValues(epochIndex, :)', chanlocs, ...
            'electrodes', 'on', ...
            'style', 'both');
        clim(lastAxes, mapLimits);
        title(sprintf('%.2f-%.2f s', ...
            epochStartTimeSec(epochIndex), epochEndTimeSec(epochIndex)));
    end

    colormap(epochFigure, jet(256));

    lastAxesPosition = lastAxes.Position;
    epochColorbar = colorbar(lastAxes, 'eastoutside');
    lastAxes.Position = lastAxesPosition;
    epochColorbar.Position = [0.935 0.15 0.018 0.70];
    epochColorbar.Label.String = colorbarLabel;
    sgtitle(epochFigure, sprintf( ...
        'Gamma-band Topography Over Time (%.1f-%.1f Hz; %.1f-s epochs)', ...
        gammaBandHz(1), gammaBandHz(2), epochLengthSec));
    exportgraphics( ...
        epochFigure, epochOutputFilePath, ...
        'Resolution', outputResolutionDpi);
end

%% Return all calculated values
GammaTopoplot_Result = struct();
GammaTopoplot_Result.sourceFilePath = CombinedPSD_Result.filePath;
GammaTopoplot_Result.wholeSpanOutputFilePath = wholeSpanOutputFilePath;
GammaTopoplot_Result.epochOutputFilePath = epochOutputFilePath;
GammaTopoplot_Result.plotMode = char(plotMode);
GammaTopoplot_Result.labels = labels;
GammaTopoplot_Result.chanlocs = chanlocs;
GammaTopoplot_Result.Fs = Fs;
GammaTopoplot_Result.gammaBandHz = gammaBandHz;
GammaTopoplot_Result.analysisDurationSec = analysisDurationSec;
GammaTopoplot_Result.epochLengthSec = epochLengthSec;
GammaTopoplot_Result.epochStartTimeSec = epochStartTimeSec;
GammaTopoplot_Result.epochEndTimeSec = epochEndTimeSec;
GammaTopoplot_Result.epochCenterTimeSec = epochCenterTimeSec;
GammaTopoplot_Result.gammaPower = gammaPower;
GammaTopoplot_Result.topoplotValues = topoplotValues;
GammaTopoplot_Result.wholeSpanGammaPower = wholeSpanGammaPower;
GammaTopoplot_Result.wholeSpanTopoplotValues = wholeSpanTopoplotValues;
GammaTopoplot_Result.displayMode = char(displayMode);
GammaTopoplot_Result.commonMapLimits = commonMapLimits;
GammaTopoplot_Result.numberOfSaturatedEpochValues = ...
    numberOfSaturatedEpochValues;
GammaTopoplot_Result.epochMapLimits = mapLimits;
GammaTopoplot_Result.wholeSpanMapLimits = wholeSpanMapLimits;
assignin('base', 'GammaTopoplot_Result', GammaTopoplot_Result);

if closeFiguresAfterSave
    if ~isempty(wholeSpanFigure)
        close(wholeSpanFigure);
    end
    if ~isempty(epochFigure)
        close(epochFigure);
    end
end

%% Local functions
function chanlocs = resolveTimeTopoplotLocations(labels, sourceFilePath, channelLocationFile)
    if strlength(string(channelLocationFile)) > 0
        allLocations = readlocs(char(channelLocationFile));
        locationIndices = matchTimeTopoplotLabels({allLocations.labels}, labels);
        chanlocs = allLocations(locationIndices);
        validateTimeTopoplotCoordinates(chanlocs);
        return;
    end

    [~, ~, sourceExtension] = fileparts(sourceFilePath);
    if strcmpi(sourceExtension, '.mat')
        sourceVariables = whos('-file', sourceFilePath);
        if any(strcmp({sourceVariables.name}, 'EEG'))
            loadedData = load(sourceFilePath, 'EEG');
            if isstruct(loadedData.EEG) && ...
                    isfield(loadedData.EEG, 'chanlocs') && ...
                    ~isempty(loadedData.EEG.chanlocs)
                try
                    locationIndices = matchTimeTopoplotLabels( ...
                        {loadedData.EEG.chanlocs.labels}, labels);
                    candidateLocations = loadedData.EEG.chanlocs(locationIndices);
                    validateTimeTopoplotCoordinates(candidateLocations);
                    chanlocs = candidateLocations;
                    return;
                catch locationError
                    warning(['The source EEG.chanlocs could not be used: %s ', ...
                        'Trying EEGLAB standard locations instead.'], ...
                        locationError.message);
                end
            end
        end
    end

    standardLocationFile = findStandardTimeTopoplotLocationFile();
    lookupLabels = normalizeTimeTopoplotLabels(labels);
    chanlocs = struct('labels', lookupLabels);
    chanlocs = pop_chanedit(chanlocs, 'lookup', standardLocationFile);
    validateTimeTopoplotCoordinates(chanlocs);

    for channelIndex = 1:numel(chanlocs)
        chanlocs(channelIndex).labels = labels{channelIndex};
    end
end

function standardLocationFile = findStandardTimeTopoplotLocationFile()
    eeglabFile = which('eeglab');
    if isempty(eeglabFile)
        error('The EEGLAB installation folder could not be located.');
    end

    eeglabFolder = fileparts(eeglabFile);
    standardLocationFile = fullfile( ...
        eeglabFolder, ...
        'functions', ...
        'supportfiles', ...
        'Standard-10-5-Cap385.sfp');

    if exist(standardLocationFile, 'file') == 2
        return;
    end

    matchingFiles = dir(fullfile( ...
        eeglabFolder, '**', 'Standard-10-5-Cap385.sfp'));
    if isempty(matchingFiles)
        error(['EEGLAB''s Standard-10-5-Cap385.sfp file was not found. ', ...
            'Set channelLocationFile to your electrode-location file.']);
    end

    standardLocationFile = fullfile( ...
        matchingFiles(1).folder, matchingFiles(1).name);
end

function indices = matchTimeTopoplotLabels(availableLabels, requestedLabels)
    availableNormalized = lower(string(normalizeTimeTopoplotLabels(availableLabels)));
    requestedNormalized = lower(string(normalizeTimeTopoplotLabels(requestedLabels)));
    indices = zeros(1, numel(requestedNormalized));

    for channelIndex = 1:numel(requestedNormalized)
        match = find(availableNormalized == requestedNormalized(channelIndex), 1);
        if isempty(match)
            error('No channel location was found for "%s".', ...
                requestedLabels{channelIndex});
        end
        indices(channelIndex) = match;
    end
end

function normalizedLabels = normalizeTimeTopoplotLabels(labels)
    normalizedLabels = cellstr(string(labels));

    for channelIndex = 1:numel(normalizedLabels)
        label = strtrim(normalizedLabels{channelIndex});
        label = regexprep(label, '^EEG[\s_-]*', '', 'ignorecase');
        label = regexprep( ...
            label, ...
            '([_-](AVG|REF|LE|A1|A2|M1|M2))$', ...
            '', ...
            'ignorecase');
        normalizedLabels{channelIndex} = label;
    end
end

function validateTimeTopoplotCoordinates(chanlocs)
    hasCoordinates = false(1, numel(chanlocs));

    for channelIndex = 1:numel(chanlocs)
        channel = chanlocs(channelIndex);
        hasPolarCoordinates = isfield(channel, 'theta') && ...
            isfield(channel, 'radius') && ...
            ~isempty(channel.theta) && ...
            ~isempty(channel.radius);
        hasCartesianCoordinates = isfield(channel, 'X') && ...
            isfield(channel, 'Y') && ...
            ~isempty(channel.X) && ...
            ~isempty(channel.Y);
        hasCoordinates(channelIndex) = ...
            hasPolarCoordinates || hasCartesianCoordinates;
    end

    if ~all(hasCoordinates)
        missingLabels = string({chanlocs(~hasCoordinates).labels});
        error(['No usable electrode coordinates were found for: %s. ', ...
            'Set channelLocationFile to the correct location file.'], ...
            strjoin(missingLabels, ', '));
    end
end

function safeName = sanitizeTimeTopoplotFileName(value)
    safeName = char(string(value));
    safeName = regexprep(safeName, '[<>:"/\\|?*]', '_');
    safeName = regexprep(safeName, '\s+', '_');
    safeName = regexprep(safeName, '_+', '_');
    safeName = strtrim(safeName);

    if isempty(safeName)
        safeName = 'EEG';
    end
end
