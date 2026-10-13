classdef PhaseConnectivityApp < matlab.apps.AppBase
    %PHASECONNECTIVITYAPP Interactive PLV/PLI/wPLI, PSD, and scalp viewer.

    properties (Access = public)
        UIFigure matlab.ui.Figure
    end

    properties (Access = private)
        FilePathField matlab.ui.control.EditField
        BrowseButton matlab.ui.control.Button
        MetadataLabel matlab.ui.control.Label
        StartField matlab.ui.control.NumericEditField
        EndField matlab.ui.control.NumericEditField
        ExportButton matlab.ui.control.Button
        StatusLabel matlab.ui.control.Label
        TabGroup matlab.ui.container.TabGroup

        ChannelTable matlab.ui.control.Table
        ConnBandEnabled matlab.ui.control.CheckBox
        ConnLowField matlab.ui.control.NumericEditField
        ConnHighField matlab.ui.control.NumericEditField
        ConnOrderField matlab.ui.control.NumericEditField
        ConnNotchDropDown matlab.ui.control.DropDown
        WindowLengthField matlab.ui.control.NumericEditField
        WindowStepField matlab.ui.control.NumericEditField
        MontageButton matlab.ui.control.Button
        MontageLabel matlab.ui.control.Label
        RunConnectivityButton matlab.ui.control.Button
        ExportWindowsCheckBox matlab.ui.control.CheckBox
        ExportSweepCheckBox matlab.ui.control.CheckBox
        MetricDropDown matlab.ui.control.DropDown
        ScopeDropDown matlab.ui.control.DropDown
        ViewDropDown matlab.ui.control.DropDown
        ThresholdLabel matlab.ui.control.Label
        ThresholdSlider matlab.ui.control.Slider
        WindowSlider matlab.ui.control.Slider
        WindowSpinner matlab.ui.control.Spinner
        WindowInfoLabel matlab.ui.control.Label
        MappingInfoLabel matlab.ui.control.Label
        ConnectivityAxes matlab.ui.control.UIAxes

        PSDChannelDropDown matlab.ui.control.DropDown
        PSDBandEnabled matlab.ui.control.CheckBox
        PSDLowField matlab.ui.control.NumericEditField
        PSDHighField matlab.ui.control.NumericEditField
        PSDOrderField matlab.ui.control.NumericEditField
        PSDNotchDropDown matlab.ui.control.DropDown
        PSDWindowField matlab.ui.control.NumericEditField
        PSDOverlapField matlab.ui.control.NumericEditField
        PSDMinFrequencyField matlab.ui.control.NumericEditField
        PSDMaxFrequencyField matlab.ui.control.NumericEditField
        PSDSourceDropDown matlab.ui.control.DropDown
        PSDScaleDropDown matlab.ui.control.DropDown
        RunPSDButton matlab.ui.control.Button
        PSDAxes matlab.ui.control.UIAxes

        PSDTopoLowField matlab.ui.control.NumericEditField
        PSDTopoHighField matlab.ui.control.NumericEditField
        PSDTopoEpochField matlab.ui.control.NumericEditField
        PSDTopoSourceDropDown matlab.ui.control.DropDown
        PSDTopoScaleDropDown matlab.ui.control.DropDown
        PSDTopoScopeDropDown matlab.ui.control.DropDown
        PSDTopoEpochSlider matlab.ui.control.Slider
        PSDTopoEpochSpinner matlab.ui.control.Spinner
        PSDTopoEpochInfoLabel matlab.ui.control.Label
        PSDTopoMappingInfoLabel matlab.ui.control.Label
        PSDTopoMontageLabel matlab.ui.control.Label
        PSDTopoExportAllCheckBox matlab.ui.control.CheckBox
        RunPSDTopoButton matlab.ui.control.Button
        PSDTopoAxes matlab.ui.control.UIAxes

        Recording = []
        ConnectivityResults = []
        PSDResults = []
        PSDTopographyResults = []
        ElectrodeMapping = []
        PSDTopographyMapping = []
        MontageSource = []
    end

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Name = 'EEG Phase Connectivity';
            app.UIFigure.Position = [80 80 1450 880];
            app.UIFigure.Color = [0.96 0.97 0.98];

            root = uigridlayout(app.UIFigure, [3 1]);
            root.RowHeight = {142, '1x', 30};
            root.Padding = [10 10 10 8];
            root.RowSpacing = 8;

            filePanel = uipanel(root, 'Title', 'Recording and shared interval');
            filePanel.Layout.Row = 1;
            fileGrid = uigridlayout(filePanel, [3 8]);
            fileGrid.RowHeight = {28, 28, 28};
            fileGrid.ColumnWidth = {72, '1x', 82, 80, 82, 80, 115, 110};
            fileGrid.Padding = [8 5 8 5];

            fieldLabel = uilabel(fileGrid, 'Text', 'EEG file');
            fieldLabel.Layout.Row = 1;
            app.FilePathField = uieditfield(fileGrid, 'text', 'Editable', 'off');
            app.FilePathField.Layout.Row = 1;
            app.FilePathField.Layout.Column = [2 6];
            app.BrowseButton = uibutton(fileGrid, 'push', 'Text', 'Browse...');
            app.BrowseButton.Layout.Row = 1;
            app.BrowseButton.Layout.Column = 7;
            app.BrowseButton.ButtonPushedFcn = @(~,~) browseFile(app);
            app.ExportButton = uibutton(fileGrid, 'push', 'Text', 'Export results');
            app.ExportButton.Layout.Row = 1;
            app.ExportButton.Layout.Column = 8;
            app.ExportButton.Enable = 'off';
            app.ExportButton.ButtonPushedFcn = @(~,~) exportResults(app);

            app.MetadataLabel = uilabel(fileGrid, 'Text', 'No recording loaded');
            app.MetadataLabel.Layout.Row = 2;
            app.MetadataLabel.Layout.Column = [1 8];
            app.MetadataLabel.FontColor = [0.25 0.30 0.36];

            fieldLabel = uilabel(fileGrid, 'Text', 'Start (s)');
            fieldLabel.Layout.Row = 3;
            app.StartField = uieditfield(fileGrid, 'numeric', 'Value', 0, ...
                'LowerLimitInclusive', 'on', 'Limits', [0 Inf]);
            app.StartField.Layout.Row = 3;
            app.StartField.Layout.Column = 2;
            fieldLabel = uilabel(fileGrid, 'Text', 'End (s)');
            fieldLabel.Layout.Row = 3;
            fieldLabel.Layout.Column = 3;
            app.EndField = uieditfield(fileGrid, 'numeric', 'Value', 10, ...
                'LowerLimitInclusive', 'off', 'Limits', [0 Inf]);
            app.EndField.Layout.Row = 3;
            app.EndField.Layout.Column = 4;
            app.StartField.ValueChangedFcn = @(~,~) invalidateShared(app);
            app.EndField.ValueChangedFcn = @(~,~) invalidateShared(app);
            helpLabel = uilabel(fileGrid, 'Text', ...
                'The same half-open interval [start, end) is used by all analyses.');
            helpLabel.Layout.Row = 3;
            helpLabel.Layout.Column = [5 8];
            helpLabel.FontAngle = 'italic';
            helpLabel.FontColor = [0.35 0.39 0.44];

            app.TabGroup = uitabgroup(root);
            app.TabGroup.Layout.Row = 2;
            connectivityTab = uitab(app.TabGroup, 'Title', 'Connectivity');
            psdTab = uitab(app.TabGroup, 'Title', 'Single-channel PSD');
            psdTopographyTab = uitab(app.TabGroup, 'Title', 'PSD Topoplot');
            createConnectivityTab(app, connectivityTab);
            createPSDTab(app, psdTab);
            createPSDTopographyTab(app, psdTopographyTab);

            app.StatusLabel = uilabel(root, 'Text', 'Ready');
            app.StatusLabel.Layout.Row = 3;
            app.StatusLabel.FontColor = [0.18 0.29 0.38];
            app.UIFigure.Visible = 'on';
        end

        function createConnectivityTab(app, parent)
            grid = uigridlayout(parent, [1 2]);
            grid.ColumnWidth = {330, '1x'};
            grid.Padding = [8 8 8 8];

            controlPanel = uipanel(grid, 'Title', 'Connectivity settings', 'Scrollable', 'on');
            controls = uigridlayout(controlPanel, [18 2]);
            controls.RowHeight = repmat({26}, 1, 18);
            controls.ColumnWidth = {'1x','1x'};
            controls.Padding = [8 8 8 8];
            controls.RowSpacing = 4;

            app.ChannelTable = uitable(controls);
            app.ChannelTable.Layout.Row = [1 5];
            app.ChannelTable.Layout.Column = [1 2];
            app.ChannelTable.ColumnName = {'Use','Channel'};
            app.ChannelTable.ColumnEditable = [true false];
            app.ChannelTable.Data = table(false(0,1), strings(0,1), ...
                'VariableNames', {'Use','Channel'});
            app.ChannelTable.CellEditCallback = @(~,~) invalidateConnectivity(app);

            app.ConnBandEnabled = uicheckbox(controls, 'Text', 'Band-pass enabled', 'Value', true);
            app.ConnBandEnabled.Layout.Row = 6;
            app.ConnBandEnabled.Layout.Column = [1 2];
            fieldLabel = uilabel(controls, 'Text', 'Low / high (Hz)');
            fieldLabel.Layout.Row = 7;
            bandGrid = uigridlayout(controls, [1 2]);
            bandGrid.Layout.Row = 7;
            bandGrid.Layout.Column = 2;
            bandGrid.Padding = [0 0 0 0];
            app.ConnLowField = uieditfield(bandGrid, 'numeric', 'Value', 25, 'Limits', [eps Inf]);
            app.ConnHighField = uieditfield(bandGrid, 'numeric', 'Value', 50, 'Limits', [eps Inf]);
            fieldLabel = uilabel(controls, 'Text', 'FIR order');
            fieldLabel.Layout.Row = 8;
            app.ConnOrderField = uieditfield(controls, 'numeric', 'Value', 20, ...
                'RoundFractionalValues', 'on', 'Limits', [1 Inf]);
            app.ConnOrderField.Layout.Row = 8;
            app.ConnOrderField.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Notch');
            fieldLabel.Layout.Row = 9;
            app.ConnNotchDropDown = uidropdown(controls, 'Items', {'Off','50 Hz','60 Hz'}, ...
                'ItemsData', [0 50 60], 'Value', 0);
            app.ConnNotchDropDown.Layout.Row = 9;
            app.ConnNotchDropDown.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Window / step (s)');
            fieldLabel.Layout.Row = 10;
            winGrid = uigridlayout(controls, [1 2]);
            winGrid.Layout.Row = 10;
            winGrid.Layout.Column = 2;
            winGrid.Padding = [0 0 0 0];
            app.WindowLengthField = uieditfield(winGrid, 'numeric', 'Value', .5, 'Limits', [eps Inf]);
            app.WindowStepField = uieditfield(winGrid, 'numeric', 'Value', .1, 'Limits', [eps Inf]);

            app.MontageButton = uibutton(controls, 'push', 'Text', 'Import montage...');
            app.MontageButton.Layout.Row = 11;
            app.MontageButton.Layout.Column = [1 2];
            app.MontageButton.ButtonPushedFcn = @(~,~) importMontage(app);
            app.MontageLabel = uilabel(controls, 'Text', 'Built-in 10-20 montage', ...
                'WordWrap', 'on');
            app.MontageLabel.Layout.Row = [12 13];
            app.MontageLabel.Layout.Column = [1 2];
            app.ExportWindowsCheckBox = uicheckbox(controls, ...
                'Text', 'Export every sliding window', 'Value', false);
            app.ExportWindowsCheckBox.Layout.Row = 14;
            app.ExportWindowsCheckBox.Layout.Column = [1 2];
            app.ExportSweepCheckBox = uicheckbox(controls, ...
                'Text', 'Export topoplot thresholds 0:0.1:0.9', 'Value', false);
            app.ExportSweepCheckBox.Layout.Row = 15;
            app.ExportSweepCheckBox.Layout.Column = [1 2];
            app.RunConnectivityButton = uibutton(controls, 'push', ...
                'Text', 'Run connectivity', 'FontWeight', 'bold', 'Enable', 'off');
            app.RunConnectivityButton.Layout.Row = 16;
            app.RunConnectivityButton.Layout.Column = [1 2];
            app.RunConnectivityButton.ButtonPushedFcn = @(~,~) runConnectivity(app);

            changed = @(~,~) invalidateConnectivity(app);
            app.ConnBandEnabled.ValueChangedFcn = changed;
            app.ConnLowField.ValueChangedFcn = changed;
            app.ConnHighField.ValueChangedFcn = changed;
            app.ConnOrderField.ValueChangedFcn = changed;
            app.ConnNotchDropDown.ValueChangedFcn = changed;
            app.WindowLengthField.ValueChangedFcn = changed;
            app.WindowStepField.ValueChangedFcn = changed;

            resultPanel = uipanel(grid, 'Title', 'Connectivity result');
            resultGrid = uigridlayout(resultPanel, [4 6]);
            resultGrid.RowHeight = {30, 42, 24, '1x'};
            resultGrid.ColumnWidth = {70, 120, 70, 90, '1x', 80};
            resultGrid.Padding = [8 6 8 8];
            fieldLabel = uilabel(resultGrid, 'Text', 'Metric');
            fieldLabel.Layout.Row = 1;
            app.MetricDropDown = uidropdown(resultGrid, 'Items', {'PLV','PLI','wPLI'}, 'Value', 'PLV');
            app.MetricDropDown.Layout.Row = 1; app.MetricDropDown.Layout.Column = 2;
            fieldLabel = uilabel(resultGrid, 'Text', 'Scope');
            fieldLabel.Layout.Row = 1;
            fieldLabel.Layout.Column = 3;
            app.ScopeDropDown = uidropdown(resultGrid, 'Items', {'Whole','Sliding'}, 'Value', 'Whole');
            app.ScopeDropDown.Layout.Row = 1; app.ScopeDropDown.Layout.Column = 4;
            app.ViewDropDown = uidropdown(resultGrid, 'Items', {'Heatmap','Scalp connections'}, 'Value', 'Heatmap');
            app.ViewDropDown.Layout.Row = 1; app.ViewDropDown.Layout.Column = [5 6];
            app.MetricDropDown.ValueChangedFcn = @(~,~) refreshConnectivityPlot(app);
            app.ScopeDropDown.ValueChangedFcn = @(~,~) refreshConnectivityPlot(app);
            app.ViewDropDown.ValueChangedFcn = @(~,~) refreshConnectivityPlot(app);

            app.ThresholdLabel = uilabel(resultGrid, 'Text', 'Threshold 0.50');
            app.ThresholdLabel.Layout.Row = 2; app.ThresholdLabel.Layout.Column = 1;
            app.ThresholdLabel.Visible = 'off';
            app.ThresholdSlider = uislider(resultGrid, 'Limits', [0 1], 'Value', .5, ...
                'MajorTicks', 0:.1:1);
            app.ThresholdSlider.Layout.Row = 2; app.ThresholdSlider.Layout.Column = [2 3];
            app.ThresholdSlider.Visible = 'off';
            app.ThresholdSlider.ValueChangingFcn = @(~,event) thresholdChanging(app, event.Value);
            app.ThresholdSlider.ValueChangedFcn = @(~,~) refreshConnectivityPlot(app);
            fieldLabel = uilabel(resultGrid, 'Text', 'Window');
            fieldLabel.Layout.Row = 2;
            fieldLabel.Layout.Column = 4;
            app.WindowSlider = uislider(resultGrid, 'Limits', [1 2], 'Value', 1, 'MajorTicks', []);
            app.WindowSlider.Layout.Row = 2; app.WindowSlider.Layout.Column = 5;
            app.WindowSpinner = uispinner(resultGrid, 'Limits', [1 1], 'Value', 1, ...
                'Step', 1, 'RoundFractionalValues', 'on');
            app.WindowSpinner.Layout.Row = 2; app.WindowSpinner.Layout.Column = 6;
            app.WindowSlider.ValueChangingFcn = @(~,event) windowChanging(app, event.Value);
            app.WindowSlider.ValueChangedFcn = @(~,~) windowChanged(app, app.WindowSlider.Value);
            app.WindowSpinner.ValueChangedFcn = @(~,~) windowChanged(app, app.WindowSpinner.Value);

            app.WindowInfoLabel = uilabel(resultGrid, 'Text', 'Run connectivity to view results.');
            app.WindowInfoLabel.Layout.Row = 3; app.WindowInfoLabel.Layout.Column = [1 3];
            app.MappingInfoLabel = uilabel(resultGrid, 'Text', '');
            app.MappingInfoLabel.Layout.Row = 3; app.MappingInfoLabel.Layout.Column = [4 6];
            app.MappingInfoLabel.HorizontalAlignment = 'right';
            app.ConnectivityAxes = uiaxes(resultGrid);
            app.ConnectivityAxes.Layout.Row = 4;
            app.ConnectivityAxes.Layout.Column = [1 6];
            title(app.ConnectivityAxes, 'No connectivity result');
        end

        function createPSDTab(app, parent)
            grid = uigridlayout(parent, [1 2]);
            grid.ColumnWidth = {330, '1x'};
            grid.Padding = [8 8 8 8];
            panel = uipanel(grid, 'Title', 'PSD settings');
            controls = uigridlayout(panel, [14 2]);
            controls.RowHeight = repmat({28}, 1, 14);
            controls.ColumnWidth = {'1x','1x'};
            controls.Padding = [8 8 8 8];
            fieldLabel = uilabel(controls, 'Text', 'Channel');
            fieldLabel.Layout.Row = 1;
            app.PSDChannelDropDown = uidropdown(controls, 'Items', {'No channels'});
            app.PSDChannelDropDown.Layout.Row = 1; app.PSDChannelDropDown.Layout.Column = 2;
            app.PSDBandEnabled = uicheckbox(controls, 'Text', 'Band-pass enabled', 'Value', false);
            app.PSDBandEnabled.Layout.Row = 2; app.PSDBandEnabled.Layout.Column = [1 2];
            fieldLabel = uilabel(controls, 'Text', 'Low / high (Hz)');
            fieldLabel.Layout.Row = 3;
            bandGrid = uigridlayout(controls, [1 2]); bandGrid.Padding = [0 0 0 0];
            bandGrid.Layout.Row = 3; bandGrid.Layout.Column = 2;
            app.PSDLowField = uieditfield(bandGrid, 'numeric', 'Value', 1, 'Limits', [eps Inf]);
            app.PSDHighField = uieditfield(bandGrid, 'numeric', 'Value', 70, 'Limits', [eps Inf]);
            fieldLabel = uilabel(controls, 'Text', 'FIR order');
            fieldLabel.Layout.Row = 4;
            app.PSDOrderField = uieditfield(controls, 'numeric', 'Value', 20, ...
                'RoundFractionalValues', 'on', 'Limits', [1 Inf]);
            app.PSDOrderField.Layout.Row = 4; app.PSDOrderField.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Notch');
            fieldLabel.Layout.Row = 5;
            app.PSDNotchDropDown = uidropdown(controls, 'Items', {'Off','50 Hz','60 Hz'}, ...
                'ItemsData', [0 50 60], 'Value', 0);
            app.PSDNotchDropDown.Layout.Row = 5; app.PSDNotchDropDown.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Welch window (s)');
            fieldLabel.Layout.Row = 6;
            app.PSDWindowField = uieditfield(controls, 'numeric', 'Value', 1, 'Limits', [eps Inf]);
            app.PSDWindowField.Layout.Row = 6; app.PSDWindowField.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Overlap fraction');
            fieldLabel.Layout.Row = 7;
            app.PSDOverlapField = uieditfield(controls, 'numeric', 'Value', .5, 'Limits', [0 .999]);
            app.PSDOverlapField.Layout.Row = 7; app.PSDOverlapField.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Display low / high (Hz)');
            fieldLabel.Layout.Row = 8;
            frequencyGrid = uigridlayout(controls, [1 2]);
            frequencyGrid.Layout.Row = 8;
            frequencyGrid.Layout.Column = 2;
            frequencyGrid.Padding = [0 0 0 0];
            app.PSDMinFrequencyField = uieditfield(frequencyGrid, 'numeric', ...
                'Value', 0, 'Limits', [0 Inf]);
            app.PSDMaxFrequencyField = uieditfield(frequencyGrid, 'numeric', ...
                'Value', 70, 'Limits', [eps Inf]);
            fieldLabel = uilabel(controls, 'Text', 'Displayed source');
            fieldLabel.Layout.Row = 9;
            app.PSDSourceDropDown = uidropdown(controls, 'Items', {'Raw','Filtered'}, 'Value', 'Raw');
            app.PSDSourceDropDown.Layout.Row = 9; app.PSDSourceDropDown.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Power scale');
            fieldLabel.Layout.Row = 10;
            app.PSDScaleDropDown = uidropdown(controls, 'Items', {'Linear','dB'}, 'Value', 'Linear');
            app.PSDScaleDropDown.Layout.Row = 10; app.PSDScaleDropDown.Layout.Column = 2;
            app.RunPSDButton = uibutton(controls, 'push', 'Text', 'Compute PSD', ...
                'FontWeight', 'bold', 'Enable', 'off');
            app.RunPSDButton.Layout.Row = 12; app.RunPSDButton.Layout.Column = [1 2];
            app.RunPSDButton.ButtonPushedFcn = @(~,~) runPSD(app);

            app.PSDChannelDropDown.ValueChangedFcn = @(~,~) invalidatePSD(app);
            filterChanged = @(~,~) invalidatePSDFilter(app);
            app.PSDBandEnabled.ValueChangedFcn = filterChanged;
            app.PSDLowField.ValueChangedFcn = filterChanged;
            app.PSDHighField.ValueChangedFcn = filterChanged;
            app.PSDOrderField.ValueChangedFcn = filterChanged;
            app.PSDNotchDropDown.ValueChangedFcn = filterChanged;
            app.PSDWindowField.ValueChangedFcn = @(~,~) invalidatePSD(app);
            app.PSDOverlapField.ValueChangedFcn = @(~,~) invalidatePSD(app);
            app.PSDMinFrequencyField.ValueChangedFcn = @(~,~) refreshPSDPlot(app);
            app.PSDMaxFrequencyField.ValueChangedFcn = @(~,~) refreshPSDPlot(app);
            app.PSDSourceDropDown.ValueChangedFcn = @(~,~) refreshPSDPlot(app);
            app.PSDScaleDropDown.ValueChangedFcn = @(~,~) refreshPSDPlot(app);

            resultPanel = uipanel(grid, 'Title', 'Welch power spectral density');
            resultGrid = uigridlayout(resultPanel, [1 1]);
            app.PSDAxes = uiaxes(resultGrid);
            title(app.PSDAxes, 'No PSD result');
            xlabel(app.PSDAxes, 'Frequency (Hz)');
        end

        function createPSDTopographyTab(app, parent)
            grid = uigridlayout(parent, [1 2]);
            grid.ColumnWidth = {330, '1x'};
            grid.Padding = [8 8 8 8];

            panel = uipanel(grid, 'Title', 'Band-power topoplot settings', ...
                'Scrollable', 'on');
            controls = uigridlayout(panel, [14 2]);
            controls.RowHeight = repmat({28}, 1, 14);
            controls.ColumnWidth = {'1x','1x'};
            controls.Padding = [8 8 8 8];
            controls.RowSpacing = 5;

            fieldLabel = uilabel(controls, 'Text', 'Band low / high (Hz)');
            fieldLabel.Layout.Row = 1;
            bandGrid = uigridlayout(controls, [1 2]);
            bandGrid.Layout.Row = 1; bandGrid.Layout.Column = 2;
            bandGrid.Padding = [0 0 0 0];
            app.PSDTopoLowField = uieditfield(bandGrid, 'numeric', ...
                'Value', 25, 'Limits', [0 Inf]);
            app.PSDTopoHighField = uieditfield(bandGrid, 'numeric', ...
                'Value', 50, 'Limits', [eps Inf]);
            fieldLabel = uilabel(controls, 'Text', 'Epoch length (s)');
            fieldLabel.Layout.Row = 2;
            app.PSDTopoEpochField = uieditfield(controls, 'numeric', ...
                'Value', 0.5, 'Limits', [eps Inf]);
            app.PSDTopoEpochField.Layout.Row = 2;
            app.PSDTopoEpochField.Layout.Column = 2;

            filterInfo = uilabel(controls, 'Text', ...
                'Uses the band-pass and notch settings from Single-channel PSD.', ...
                'WordWrap', 'on', 'FontAngle', 'italic');
            filterInfo.Layout.Row = [3 4];
            filterInfo.Layout.Column = [1 2];
            filterInfo.FontColor = [0.35 0.39 0.44];
            montageButton = uibutton(controls, 'push', 'Text', 'Import montage...');
            montageButton.Layout.Row = 5;
            montageButton.Layout.Column = [1 2];
            montageButton.ButtonPushedFcn = @(~,~) importMontage(app);
            app.PSDTopoMontageLabel = uilabel(controls, ...
                'Text', 'Built-in 10-20 montage', 'WordWrap', 'on');
            app.PSDTopoMontageLabel.Layout.Row = [6 7];
            app.PSDTopoMontageLabel.Layout.Column = [1 2];

            fieldLabel = uilabel(controls, 'Text', 'Displayed source');
            fieldLabel.Layout.Row = 8;
            app.PSDTopoSourceDropDown = uidropdown(controls, ...
                'Items', {'Raw','Filtered'}, 'Value', 'Raw');
            app.PSDTopoSourceDropDown.Layout.Row = 8;
            app.PSDTopoSourceDropDown.Layout.Column = 2;
            fieldLabel = uilabel(controls, 'Text', 'Power scale');
            fieldLabel.Layout.Row = 9;
            app.PSDTopoScaleDropDown = uidropdown(controls, ...
                'Items', {'Linear','dB'}, 'Value', 'Linear');
            app.PSDTopoScaleDropDown.Layout.Row = 9;
            app.PSDTopoScaleDropDown.Layout.Column = 2;
            app.PSDTopoExportAllCheckBox = uicheckbox(controls, ...
                'Text', 'Export tiled map of all epochs', 'Value', false);
            app.PSDTopoExportAllCheckBox.Layout.Row = 10;
            app.PSDTopoExportAllCheckBox.Layout.Column = [1 2];
            app.RunPSDTopoButton = uibutton(controls, 'push', ...
                'Text', 'Compute PSD Topoplot', 'FontWeight', 'bold', 'Enable', 'off');
            app.RunPSDTopoButton.Layout.Row = 12;
            app.RunPSDTopoButton.Layout.Column = [1 2];
            app.RunPSDTopoButton.ButtonPushedFcn = @(~,~) runPSDTopography(app);

            app.PSDTopoLowField.ValueChangedFcn = @(~,~) invalidatePSDTopography(app);
            app.PSDTopoHighField.ValueChangedFcn = @(~,~) invalidatePSDTopography(app);
            app.PSDTopoEpochField.ValueChangedFcn = @(~,~) invalidatePSDTopography(app);
            app.PSDTopoSourceDropDown.ValueChangedFcn = @(~,~) refreshPSDTopographyPlot(app);
            app.PSDTopoScaleDropDown.ValueChangedFcn = @(~,~) refreshPSDTopographyPlot(app);

            resultPanel = uipanel(grid, 'Title', 'Multichannel band-power scalp map');
            resultGrid = uigridlayout(resultPanel, [4 6]);
            resultGrid.RowHeight = {30, 42, 24, '1x'};
            resultGrid.ColumnWidth = {70, 120, 70, 90, '1x', 80};
            resultGrid.Padding = [8 6 8 8];
            fieldLabel = uilabel(resultGrid, 'Text', 'Scope');
            fieldLabel.Layout.Row = 1;
            app.PSDTopoScopeDropDown = uidropdown(resultGrid, ...
                'Items', {'Whole','Epoch'}, 'Value', 'Whole');
            app.PSDTopoScopeDropDown.Layout.Row = 1;
            app.PSDTopoScopeDropDown.Layout.Column = 2;
            app.PSDTopoScopeDropDown.ValueChangedFcn = @(~,~) refreshPSDTopographyPlot(app);
            fieldLabel = uilabel(resultGrid, 'Text', 'Epoch');
            fieldLabel.Layout.Row = 2;
            fieldLabel.Layout.Column = 1;
            app.PSDTopoEpochSlider = uislider(resultGrid, ...
                'Limits', [1 2], 'Value', 1, 'MajorTicks', [], 'Enable', 'off');
            app.PSDTopoEpochSlider.Layout.Row = 2;
            app.PSDTopoEpochSlider.Layout.Column = [2 5];
            app.PSDTopoEpochSpinner = uispinner(resultGrid, ...
                'Limits', [1 1], 'Value', 1, 'Step', 1, ...
                'RoundFractionalValues', 'on', 'Enable', 'off');
            app.PSDTopoEpochSpinner.Layout.Row = 2;
            app.PSDTopoEpochSpinner.Layout.Column = 6;
            app.PSDTopoEpochSlider.ValueChangingFcn = ...
                @(~,event) psdTopoEpochChanged(app, event.Value);
            app.PSDTopoEpochSlider.ValueChangedFcn = ...
                @(~,~) psdTopoEpochChanged(app, app.PSDTopoEpochSlider.Value);
            app.PSDTopoEpochSpinner.ValueChangedFcn = ...
                @(~,~) psdTopoEpochChanged(app, app.PSDTopoEpochSpinner.Value);
            app.PSDTopoEpochInfoLabel = uilabel(resultGrid, ...
                'Text', 'Compute a topoplot to view results.');
            app.PSDTopoEpochInfoLabel.Layout.Row = 3;
            app.PSDTopoEpochInfoLabel.Layout.Column = [1 3];
            app.PSDTopoMappingInfoLabel = uilabel(resultGrid, 'Text', '');
            app.PSDTopoMappingInfoLabel.Layout.Row = 3;
            app.PSDTopoMappingInfoLabel.Layout.Column = [4 6];
            app.PSDTopoMappingInfoLabel.HorizontalAlignment = 'right';
            app.PSDTopoAxes = uiaxes(resultGrid);
            app.PSDTopoAxes.Layout.Row = 4;
            app.PSDTopoAxes.Layout.Column = [1 6];
            title(app.PSDTopoAxes, 'No PSD topoplot result');
        end

        function browseFile(app)
            [file, folder] = uigetfile({'*.edf;*.mat','EEG files (*.edf, *.mat)'}, ...
                'Select an EEG recording');
            if isequal(file, 0)
                return;
            end
            loadRecording(app, fullfile(folder, file));
        end

        function loadRecording(app, path)
            dialog = uiprogressdlg(app.UIFigure, 'Title', 'Loading recording', ...
                'Message', 'Reading EEG channels...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() delete(dialog)); %#ok<NASGU>
            try
                recording = load_phase_connectivity_recording(path);
                app.Recording = recording;
                app.FilePathField.Value = path;
                app.StartField.Value = 0;
                app.EndField.Value = min(10, recording.duration_sec);
                app.MetadataLabel.Text = sprintf( ...
                    'Fs %.6g Hz  |  %.3f s  |  %d samples  |  %d EEG channels', ...
                    recording.Fs, recording.duration_sec, recording.num_samples, ...
                    numel(recording.channel_names));
                app.ChannelTable.Data = table(true(numel(recording.channel_names),1), ...
                    string(recording.channel_names), 'VariableNames', {'Use','Channel'});
                app.PSDChannelDropDown.Items = recording.channel_names;
                app.PSDChannelDropDown.Value = recording.channel_names{1};
                displayHigh = min(70, recording.Fs / 2);
                app.PSDMinFrequencyField.Limits = [0 max(eps, displayHigh)];
                app.PSDMinFrequencyField.Value = 0;
                app.PSDMaxFrequencyField.Value = displayHigh;
                app.PSDMaxFrequencyField.Limits = [eps recording.Fs / 2];
                topoHigh = min(50, recording.Fs / 2);
                topoLow = min(25, topoHigh / 2);
                app.PSDTopoLowField.Limits = [0 max(eps, topoHigh)];
                app.PSDTopoLowField.Value = topoLow;
                app.PSDTopoHighField.Limits = [eps recording.Fs / 2];
                app.PSDTopoHighField.Value = topoHigh;
                if ~isempty(recording.chanlocs)
                    app.MontageSource = recording.chanlocs;
                    app.MontageLabel.Text = 'Channel locations embedded in MAT recording';
                    app.PSDTopoMontageLabel.Text = ...
                        'Channel locations embedded in MAT recording';
                else
                    app.MontageSource = [];
                    app.MontageLabel.Text = 'Built-in 10-20 montage';
                    app.PSDTopoMontageLabel.Text = 'Built-in 10-20 montage';
                end
                app.RunConnectivityButton.Enable = 'on';
                app.RunPSDButton.Enable = 'on';
                app.RunPSDTopoButton.Enable = 'on';
                invalidateShared(app);
                app.StatusLabel.Text = sprintf('Loaded %s', recording.file_name);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'Unable to load recording');
                app.StatusLabel.Text = 'Load failed';
            end
        end

        function invalidateShared(app)
            invalidateConnectivity(app);
            invalidatePSD(app);
            invalidatePSDTopography(app);
        end

        function invalidateConnectivity(app)
            app.ConnectivityResults = [];
            app.ElectrodeMapping = [];
            cla(app.ConnectivityAxes);
            title(app.ConnectivityAxes, 'Connectivity settings changed — run again');
            updateExportState(app);
        end

        function invalidatePSD(app)
            app.PSDResults = [];
            cla(app.PSDAxes);
            title(app.PSDAxes, 'PSD settings changed — compute again');
            updateExportState(app);
        end

        function invalidatePSDFilter(app)
            invalidatePSD(app);
            invalidatePSDTopography(app);
        end

        function invalidatePSDTopography(app)
            app.PSDTopographyResults = [];
            app.PSDTopographyMapping = [];
            cla(app.PSDTopoAxes);
            title(app.PSDTopoAxes, 'PSD topoplot settings changed — compute again');
            app.PSDTopoEpochInfoLabel.Text = 'Compute a topoplot to view results.';
            app.PSDTopoMappingInfoLabel.Text = '';
            app.PSDTopoEpochSlider.Enable = 'off';
            app.PSDTopoEpochSpinner.Enable = 'off';
            updateExportState(app);
        end

        function settings = connectivityFilterSettings(app)
            settings = struct('bandpass_enabled', app.ConnBandEnabled.Value, ...
                'band_hz', [app.ConnLowField.Value app.ConnHighField.Value], ...
                'bandpass_order', app.ConnOrderField.Value, ...
                'notch_hz', app.ConnNotchDropDown.Value);
        end

        function settings = psdFilterSettings(app)
            settings = struct('bandpass_enabled', app.PSDBandEnabled.Value, ...
                'band_hz', [app.PSDLowField.Value app.PSDHighField.Value], ...
                'bandpass_order', app.PSDOrderField.Value, ...
                'notch_hz', app.PSDNotchDropDown.Value);
        end

        function runConnectivity(app)
            if isempty(app.Recording), return; end
            selected = app.ChannelTable.Data.Channel(app.ChannelTable.Data.Use);
            if numel(selected) < 2
                uialert(app.UIFigure, 'Select at least two channels.', 'Channel selection');
                return;
            end
            estimated = estimateConnectivityBytes(app, numel(selected));
            if estimated > 1024^3
                answer = uiconfirm(app.UIFigure, sprintf( ...
                    'Sliding matrices may require approximately %.2f GiB. Continue?', ...
                    estimated / 1024^3), 'Large analysis', ...
                    'Options', {'Continue','Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
                if ~strcmp(answer, 'Continue'), return; end
            end
            dialog = uiprogressdlg(app.UIFigure, 'Title', 'Connectivity analysis', ...
                'Message', 'Filtering and calculating PLV, PLI, and wPLI...', ...
                'Indeterminate', 'on');
            cleanup = onCleanup(@() delete(dialog)); %#ok<NASGU>
            try
                options = struct('start_time_sec', app.StartField.Value, ...
                    'end_time_sec', app.EndField.Value, ...
                    'channel_names', {cellstr(selected)}, ...
                    'window_length_sec', app.WindowLengthField.Value, ...
                    'window_step_sec', app.WindowStepField.Value, ...
                    'filter_settings', connectivityFilterSettings(app));
                app.ConnectivityResults = run_phase_connectivity_analysis(app.Recording, options);
                app.ElectrodeMapping = resolve_electrode_locations( ...
                    app.ConnectivityResults.channel_names, app.MontageSource);
                configureWindowControls(app);
                refreshConnectivityPlot(app);
                app.StatusLabel.Text = sprintf('Connectivity complete: %d windows', ...
                    numel(app.ConnectivityResults.time_axis));
                updateExportState(app);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'Connectivity analysis failed');
                app.StatusLabel.Text = 'Connectivity analysis failed';
            end
        end

        function bytes = estimateConnectivityBytes(app, numChannels)
            Fs = app.Recording.Fs;
            samples = floor(app.EndField.Value * Fs) - floor(app.StartField.Value * Fs);
            window = round(app.WindowLengthField.Value * Fs);
            step = max(1, round(app.WindowStepField.Value * Fs));
            windows = max(0, floor((samples - window) / step) + 1);
            bytes = 3 * numChannels^2 * windows * 8;
        end

        function configureWindowControls(app)
            count = numel(app.ConnectivityResults.time_axis);
            app.WindowSpinner.Limits = [1 max(1, count)];
            app.WindowSpinner.Value = 1;
            app.WindowSlider.Limits = [1 max(2, count)];
            app.WindowSlider.Value = 1;
            if count <= 1
                app.WindowSlider.Enable = 'off';
                app.WindowSpinner.Enable = 'off';
            else
                app.WindowSlider.Enable = 'on';
                app.WindowSpinner.Enable = 'on';
            end
        end

        function [matrix, description] = currentConnectivityMatrix(app)
            results = app.ConnectivityResults;
            metric = app.MetricDropDown.Value;
            prefix = lower(metric);
            if strcmp(metric, 'wPLI'), prefix = 'wpli'; end
            if strcmp(app.ScopeDropDown.Value, 'Whole')
                matrix = results.([prefix '_matrix']);
                description = sprintf('Whole segment %.3f-%.3f s', ...
                    results.start_time_sec, results.end_time_sec);
            else
                index = min(numel(results.time_axis), max(1, round(app.WindowSpinner.Value)));
                values = results.([prefix '_3d_matrix']);
                matrix = values(:, :, index);
                description = sprintf('Window %d | center %.4f s', index, ...
                    results.absolute_time_axis(index));
            end
        end

        function refreshConnectivityPlot(app)
            if isempty(app.ConnectivityResults), return; end
            [matrix, description] = currentConnectivityMatrix(app);
            metric = app.MetricDropDown.Value;
            app.ThresholdLabel.Text = sprintf('Threshold %.2f', app.ThresholdSlider.Value);
            sliding = strcmp(app.ScopeDropDown.Value, 'Sliding');
            app.WindowSlider.Enable = matlab.lang.OnOffSwitchState(sliding && ...
                numel(app.ConnectivityResults.time_axis) > 1);
            app.WindowSpinner.Enable = app.WindowSlider.Enable;
            if sliding
                app.WindowInfoLabel.Text = description;
            else
                app.WindowInfoLabel.Text = 'Whole-segment result';
            end

            if strcmp(app.ViewDropDown.Value, 'Heatmap')
                app.ThresholdLabel.Visible = 'off';
                app.ThresholdSlider.Visible = 'off';
                app.MappingInfoLabel.Text = '';
                displayed = matrix;
                displayed(triu(true(size(displayed)), 1)) = NaN;
                cla(app.ConnectivityAxes);
                imageHandle = imagesc(app.ConnectivityAxes, displayed, [0 1]);
                imageHandle.AlphaData = ~isnan(displayed);
                app.ConnectivityAxes.Color = [1 1 1];
                app.ConnectivityAxes.YDir = 'normal';
                axis(app.ConnectivityAxes, 'square');
                names = app.ConnectivityResults.channel_names;
                app.ConnectivityAxes.XTick = 1:numel(names);
                app.ConnectivityAxes.YTick = 1:numel(names);
                app.ConnectivityAxes.XTickLabel = names;
                app.ConnectivityAxes.YTickLabel = names;
                app.ConnectivityAxes.XTickLabelRotation = 45;
                app.ConnectivityAxes.FontSize = 8;
                colormap(app.ConnectivityAxes, jet(256));
                colorbar(app.ConnectivityAxes);
                title(app.ConnectivityAxes, sprintf('%s | %s', metric, description), ...
                    'Interpreter', 'none');
            else
                app.ThresholdLabel.Visible = 'on';
                app.ThresholdSlider.Visible = 'on';
                mapping = app.ElectrodeMapping;
                if isempty(mapping) || mapping.num_mapped < 2
                    cla(app.ConnectivityAxes);
                    title(app.ConnectivityAxes, 'At least two mapped channels are required.');
                    return;
                end
                mapped = matrix(mapping.analysis_indices, mapping.analysis_indices);
                threshold = app.ThresholdSlider.Value;
                mask = tril(true(size(mapped)), -1) & mapped >= threshold & mapped > 0;
                [first, second] = find(mask);
                displayStruct = struct('chanPairs', [first second], ...
                    'connectStrength', mapped(mask), 'connectStrengthLimits', [0 1]);
                stats = topoplot_connect(displayStruct, mapping, app.ConnectivityAxes);
                title(app.ConnectivityAxes, sprintf('%s | %s | threshold %.2f', ...
                    metric, description, threshold), 'Interpreter', 'none');
                app.MappingInfoLabel.Text = sprintf('%d/%d mapped | %d connections', ...
                    mapping.num_mapped, numel(app.ConnectivityResults.channel_names), ...
                    stats.numConnections);
            end
        end

        function thresholdChanging(app, value)
            app.ThresholdSlider.Value = value;
            refreshConnectivityPlot(app);
        end

        function windowChanging(app, value)
            windowChanged(app, round(value));
        end

        function windowChanged(app, value)
            if isempty(app.ConnectivityResults), return; end
            count = numel(app.ConnectivityResults.time_axis);
            value = min(count, max(1, round(value)));
            app.WindowSpinner.Value = value;
            app.WindowSlider.Value = value;
            refreshConnectivityPlot(app);
        end

        function importMontage(app)
            if isempty(app.Recording), return; end
            [file, folder] = uigetfile({'*.set;*.mat','Channel locations (*.set, *.mat)'}, ...
                'Select a montage source');
            if isequal(file, 0), return; end
            source = fullfile(folder, file);
            try
                recordingMapping = resolve_electrode_locations( ...
                    app.Recording.channel_names, source);
                app.MontageSource = source;
                if ~isempty(app.ConnectivityResults)
                    app.ElectrodeMapping = resolve_electrode_locations( ...
                        app.ConnectivityResults.channel_names, source);
                end
                if ~isempty(app.PSDTopographyResults)
                    app.PSDTopographyMapping = resolve_electrode_locations( ...
                        app.PSDTopographyResults.channel_names, source);
                else
                    app.PSDTopographyMapping = recordingMapping;
                end
                labelText = sprintf('%s | %d/%d mapped', file, ...
                    recordingMapping.num_mapped, numel(app.Recording.channel_names));
                app.MontageLabel.Text = labelText;
                app.PSDTopoMontageLabel.Text = labelText;
                refreshConnectivityPlot(app);
                refreshPSDTopographyPlot(app);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'Montage import failed');
            end
        end

        function runPSD(app)
            if isempty(app.Recording), return; end
            dialog = uiprogressdlg(app.UIFigure, 'Title', 'PSD analysis', ...
                'Message', 'Calculating raw and filtered Welch PSD...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() delete(dialog)); %#ok<NASGU>
            try
                options = struct('start_time_sec', app.StartField.Value, ...
                    'end_time_sec', app.EndField.Value, ...
                    'channel_name', app.PSDChannelDropDown.Value, ...
                    'window_length_sec', app.PSDWindowField.Value, ...
                    'overlap_fraction', app.PSDOverlapField.Value, ...
                    'max_frequency_hz', app.Recording.Fs / 2, ...
                    'filter_settings', psdFilterSettings(app));
                app.PSDResults = compute_single_channel_psd(app.Recording, options);
                refreshPSDPlot(app);
                app.StatusLabel.Text = sprintf('PSD complete: %s', app.PSDResults.channel_name);
                updateExportState(app);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'PSD analysis failed');
                app.StatusLabel.Text = 'PSD analysis failed';
            end
        end

        function refreshPSDPlot(app)
            if isempty(app.PSDResults), return; end
            lowFrequency = app.PSDMinFrequencyField.Value;
            highFrequency = app.PSDMaxFrequencyField.Value;
            nyquist = app.PSDResults.Fs / 2;
            if lowFrequency < 0 || highFrequency <= lowFrequency || highFrequency > nyquist
                cla(app.PSDAxes);
                title(app.PSDAxes, sprintf( ...
                    'Display range must satisfy 0 <= low < high <= %.3g Hz.', nyquist));
                return;
            end
            source = app.PSDSourceDropDown.Value;
            if strcmp(source, 'Filtered')
                values = app.PSDResults.filtered_psd;
            else
                values = app.PSDResults.raw_psd;
            end
            if strcmp(app.PSDScaleDropDown.Value, 'dB')
                values = 10 * log10(values + eps);
                unit = 'PSD (dB/Hz)';
            else
                unit = 'PSD (\muV^2/Hz)';
            end
            visible = app.PSDResults.frequency >= lowFrequency & ...
                app.PSDResults.frequency <= highFrequency;
            if nnz(visible) < 2
                cla(app.PSDAxes);
                title(app.PSDAxes, 'The selected display range contains too few PSD bins.');
                return;
            end
            cla(app.PSDAxes);
            plot(app.PSDAxes, app.PSDResults.frequency(visible), values(visible), ...
                'LineWidth', 1.35, 'Color', [0.04 0.36 0.58]);
            grid(app.PSDAxes, 'on');
            xlim(app.PSDAxes, [lowFrequency highFrequency]);
            app.PSDAxes.YLimMode = 'auto';
            xlabel(app.PSDAxes, 'Frequency (Hz)');
            ylabel(app.PSDAxes, unit);
            title(app.PSDAxes, sprintf('%s | %s PSD | %.3f-%.3f s', ...
                app.PSDResults.channel_name, source, app.PSDResults.start_time_sec, ...
                app.PSDResults.end_time_sec), 'Interpreter', 'none');
        end

        function runPSDTopography(app)
            if isempty(app.Recording), return; end
            dialog = uiprogressdlg(app.UIFigure, 'Title', 'PSD topography', ...
                'Message', 'Calculating multichannel band power...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() delete(dialog)); %#ok<NASGU>
            try
                options = struct('start_time_sec', app.StartField.Value, ...
                    'end_time_sec', app.EndField.Value, ...
                    'band_hz', [app.PSDTopoLowField.Value app.PSDTopoHighField.Value], ...
                    'epoch_length_sec', app.PSDTopoEpochField.Value, ...
                    'filter_settings', psdFilterSettings(app));
                app.PSDTopographyResults = compute_psd_topography(app.Recording, options);
                app.PSDTopographyMapping = resolve_electrode_locations( ...
                    app.PSDTopographyResults.channel_names, app.MontageSource);
                configurePSDTopographyEpochControls(app);
                refreshPSDTopographyPlot(app);
                app.StatusLabel.Text = sprintf( ...
                    'PSD topoplot complete: %d epochs, %.3g-%.3g Hz', ...
                    app.PSDTopographyResults.num_epochs, ...
                    app.PSDTopographyResults.band_hz(1), ...
                    app.PSDTopographyResults.band_hz(2));
                updateExportState(app);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'PSD topography failed');
                app.StatusLabel.Text = 'PSD topography failed';
            end
        end

        function configurePSDTopographyEpochControls(app)
            count = app.PSDTopographyResults.num_epochs;
            app.PSDTopoEpochSpinner.Limits = [1 max(1, count)];
            app.PSDTopoEpochSpinner.Value = 1;
            app.PSDTopoEpochSlider.Limits = [1 max(2, count)];
            app.PSDTopoEpochSlider.Value = 1;
        end

        function psdTopoEpochChanged(app, value)
            if isempty(app.PSDTopographyResults), return; end
            count = app.PSDTopographyResults.num_epochs;
            value = min(count, max(1, round(value)));
            app.PSDTopoEpochSpinner.Value = value;
            app.PSDTopoEpochSlider.Value = value;
            refreshPSDTopographyPlot(app);
        end

        function [values, allValues, description, colorbarLabel] = ...
                currentPSDTopographyDisplay(app)
            result = app.PSDTopographyResults;
            source = app.PSDTopoSourceDropDown.Value;
            if strcmp(source, 'Filtered')
                epochs = result.filtered_band_power;
                whole = result.filtered_whole_band_power;
            else
                epochs = result.raw_band_power;
                whole = result.raw_whole_band_power;
            end
            if strcmp(app.PSDTopoScaleDropDown.Value, 'dB')
                epochs = 10 * log10(max(epochs, realmin('double')));
                whole = 10 * log10(max(whole, realmin('double')));
                colorbarLabel = 'Integrated band power (dB)';
            else
                colorbarLabel = 'Integrated band power (signal-unit^2)';
            end
            allValues = [whole; epochs];
            if strcmp(app.PSDTopoScopeDropDown.Value, 'Whole')
                values = whole;
                description = sprintf('Whole %.3f-%.3f s', ...
                    result.start_time_sec, result.end_time_sec);
            else
                index = min(result.num_epochs, max(1, ...
                    round(app.PSDTopoEpochSpinner.Value)));
                values = epochs(index, :);
                description = sprintf('Epoch %d | %.3f-%.3f s', index, ...
                    result.epoch_start_time_sec(index), ...
                    result.epoch_end_time_sec(index));
            end
        end

        function refreshPSDTopographyPlot(app)
            if isempty(app.PSDTopographyResults), return; end
            result = app.PSDTopographyResults;
            epochView = strcmp(app.PSDTopoScopeDropDown.Value, 'Epoch');
            enabled = matlab.lang.OnOffSwitchState(epochView && result.num_epochs > 1);
            app.PSDTopoEpochSlider.Enable = enabled;
            app.PSDTopoEpochSpinner.Enable = enabled;
            [values, allValues, description, colorbarLabel] = ...
                currentPSDTopographyDisplay(app);
            mapping = app.PSDTopographyMapping;
            if isempty(mapping) || mapping.num_mapped < 3
                cla(app.PSDTopoAxes);
                title(app.PSDTopoAxes, 'At least three mapped channels are required.');
                mappedCount = 0;
                if ~isempty(mapping), mappedCount = mapping.num_mapped; end
                app.PSDTopoMappingInfoLabel.Text = sprintf('%d/%d mapped', ...
                    mappedCount, numel(result.channel_names));
                return;
            end
            displayedAll = allValues(:, mapping.analysis_indices);
            mapLimits = [min(displayedAll(:)), max(displayedAll(:))];
            try
                plot_psd_topography(values, mapping, app.PSDTopoAxes, struct( ...
                    'map_limits', mapLimits, 'colorbar_label', colorbarLabel));
                title(app.PSDTopoAxes, sprintf( ...
                    '%s | %.3g-%.3g Hz | %s | %s', ...
                    app.PSDTopoSourceDropDown.Value, result.band_hz(1), ...
                    result.band_hz(2), app.PSDTopoScaleDropDown.Value, description), ...
                    'Interpreter', 'none');
            catch errorInfo
                cla(app.PSDTopoAxes);
                title(app.PSDTopoAxes, errorInfo.message, 'Interpreter', 'none');
            end
            tailText = '';
            if result.unused_tail_sec > 0
                tailText = sprintf(' | %.4g s tail unused', result.unused_tail_sec);
            end
            app.PSDTopoEpochInfoLabel.Text = [description tailText];
            unmatched = numel(mapping.unmatched_channels);
            if unmatched > 0
                unmatchedText = strjoin(mapping.unmatched_channels, ', ');
                app.PSDTopoMappingInfoLabel.Text = sprintf( ...
                    '%d/%d mapped | unmatched: %s', mapping.num_mapped, ...
                    numel(result.channel_names), unmatchedText);
                app.PSDTopoMappingInfoLabel.Tooltip = ...
                    ['Unmatched channels: ' unmatchedText];
            else
                app.PSDTopoMappingInfoLabel.Text = sprintf('%d/%d mapped', ...
                    mapping.num_mapped, numel(result.channel_names));
                app.PSDTopoMappingInfoLabel.Tooltip = '';
            end
        end

        function updateExportState(app)
            if isempty(app.ConnectivityResults) && isempty(app.PSDResults) && ...
                    isempty(app.PSDTopographyResults)
                app.ExportButton.Enable = 'off';
            else
                app.ExportButton.Enable = 'on';
            end
        end

        function exportResults(app)
            folder = uigetdir(pwd, 'Select result output folder');
            if isequal(folder, 0), return; end
            bundle = struct('connectivity', app.ConnectivityResults, ...
                'psd', app.PSDResults, 'electrode_mapping', app.ElectrodeMapping, ...
                'psd_topography', app.PSDTopographyResults, ...
                'psd_topography_mapping', app.PSDTopographyMapping);
            options = struct('export_connectivity', ~isempty(app.ConnectivityResults), ...
                'export_all_windows', app.ExportWindowsCheckBox.Value, ...
                'export_psd', ~isempty(app.PSDResults), ...
                'psd_source', app.PSDSourceDropDown.Value, ...
                'psd_in_db', strcmp(app.PSDScaleDropDown.Value, 'dB'), ...
                'psd_frequency_range_hz', [app.PSDMinFrequencyField.Value ...
                    app.PSDMaxFrequencyField.Value], ...
                'export_topoplot', ~isempty(app.ConnectivityResults), ...
                'export_threshold_sweep', app.ExportSweepCheckBox.Value, ...
                'topoplot_metric', app.MetricDropDown.Value, ...
                'topoplot_scope', app.ScopeDropDown.Value, ...
                'topoplot_window_index', app.WindowSpinner.Value, ...
                'topoplot_threshold', app.ThresholdSlider.Value, ...
                'export_psd_topography', ~isempty(app.PSDTopographyResults), ...
                'export_all_psd_topography_epochs', app.PSDTopoExportAllCheckBox.Value, ...
                'psd_topography_source', app.PSDTopoSourceDropDown.Value, ...
                'psd_topography_in_db', strcmp(app.PSDTopoScaleDropDown.Value, 'dB'), ...
                'psd_topography_scope', app.PSDTopoScopeDropDown.Value, ...
                'psd_topography_epoch_index', app.PSDTopoEpochSpinner.Value, ...
                'resolution_dpi', 300);
            dialog = uiprogressdlg(app.UIFigure, 'Title', 'Exporting', ...
                'Message', 'Writing MAT and figure files...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() delete(dialog)); %#ok<NASGU>
            try
                files = export_phase_connectivity_results(bundle, folder, options);
                app.StatusLabel.Text = sprintf('Exported %d files to %s', numel(files), folder);
            catch errorInfo
                uialert(app.UIFigure, errorInfo.message, 'Export failed');
                app.StatusLabel.Text = 'Export failed';
            end
        end
    end

    methods (Access = public)
        function app = PhaseConnectivityApp
            createComponents(app);
            registerApp(app, app.UIFigure);
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end
end
