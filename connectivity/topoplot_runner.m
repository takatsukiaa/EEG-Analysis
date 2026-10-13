%% Automatically save PLV topoplots across thresholds
% Required variables in the workspace:
%   plv_matrix  - square channel-by-channel PLV/connectivity matrix
%   EEG         - EEGLAB EEG struct with EEG.chanlocs

clear ds chanlocs adjMat;

thresholds = 0:0.1:0.9;
outputFolder = fullfile(pwd, 'PLV_topoplot');
outputPrefix = 'No2_segment2_during';
outputResolutionDpi = 300;
showFigures = false;
closeFiguresAfterSave = true;

if ~exist('plv_matrix', 'var')
    error('plv_matrix was not found in the workspace.');
end

if ~exist('EEG', 'var') || ~isfield(EEG, 'chanlocs')
    error('EEG.chanlocs was not found in the workspace.');
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

figureVisibilityOptions = {'off', 'on'};
figureVisibility = figureVisibilityOptions{double(showFigures) + 1};

chanlocs = EEG.chanlocs;
Topoplot_Result = struct();

for thresholdIndex = 1:numel(thresholds)
    threshold = thresholds(thresholdIndex);

    lowerTriangleMask = tril(true(size(plv_matrix)), -1);
    adjMat = plv_matrix >= threshold & lowerTriangleMask;
    [chan1, chan2] = ind2sub(size(adjMat), find(adjMat));
    connectStrength = plv_matrix(adjMat);

    fig = figure('Color', 'w', 'Visible', figureVisibility);
    colormap(fig, 'jet');
    clim([0, 1]);
    colorbar;

    ds = struct();
    ds.chanPairs = [chan1, chan2];
    ds.connectStrength = connectStrength;
    ds.connectStrengthLimits = [0, 1];

    topoplot_connect(ds, chanlocs);

    if isempty(connectStrength)
        titleText = sprintf('PLV connection strengths (>=%.2f)\nNo connections', threshold);
        minimumStrength = NaN;
        maximumStrength = NaN;
    else
        minimumStrength = min(connectStrength);
        maximumStrength = max(connectStrength);
        titleText = sprintf('PLV connection strengths (>=%.2f)\nMin= %.2f, Max= %.2f', ...
            threshold, minimumStrength, maximumStrength);
    end

    title(titleText);

    outputFileName = sprintf('%s_%.1f.png', outputPrefix, threshold);
    outputFilePath = fullfile(outputFolder, outputFileName);
    exportgraphics(fig, outputFilePath, 'Resolution', outputResolutionDpi);

    Topoplot_Result(thresholdIndex).threshold = threshold;
    Topoplot_Result(thresholdIndex).outputFilePath = outputFilePath;
    Topoplot_Result(thresholdIndex).numConnections = numel(connectStrength);
    Topoplot_Result(thresholdIndex).minimumStrength = minimumStrength;
    Topoplot_Result(thresholdIndex).maximumStrength = maximumStrength;

    if closeFiguresAfterSave
        close(fig);
    end
end

assignin('base', 'Topoplot_Result', Topoplot_Result);
