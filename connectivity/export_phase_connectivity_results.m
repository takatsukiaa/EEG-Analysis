function exportedFiles = export_phase_connectivity_results(bundle, outputFolder, options)
%EXPORT_PHASE_CONNECTIVITY_RESULTS Save numerical results and 300-DPI plots.
%   BUNDLE may contain connectivity, psd, electrode_mapping, and display.

    if nargin < 3 || isempty(options)
        options = struct();
    end
    options = normalize_options(options);
    outputFolder = char(string(outputFolder));
    if isempty(outputFolder)
        error('export_phase_connectivity_results:InvalidFolder', ...
            'Select an output folder.');
    end
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    if ~isstruct(bundle) || ~isscalar(bundle)
        error('export_phase_connectivity_results:InvalidBundle', ...
            'bundle must be a scalar structure.');
    end

    matFile = fullfile(outputFolder, 'PhaseConnectivityResults.mat');
    save(matFile, 'bundle', '-v7.3');
    exportedFiles = {matFile};

    if options.export_connectivity && isfield(bundle, 'connectivity') && ...
            ~isempty(bundle.connectivity)
        connectivity = bundle.connectivity;
        wholeFile = fullfile(outputFolder, 'PhaseConnectivity_Whole.png');
        export_comparison(connectivity.plv_matrix, connectivity.pli_matrix, ...
            connectivity.wpli_matrix, connectivity.channel_names, ...
            sprintf('Whole segment %.3f-%.3f s', connectivity.start_time_sec, ...
            connectivity.end_time_sec), wholeFile, options.resolution_dpi);
        exportedFiles{end + 1} = wholeFile; %#ok<AGROW>

        if options.export_all_windows
            for index = 1:numel(connectivity.time_axis)
                outputFile = fullfile(outputFolder, sprintf( ...
                    'PhaseConnectivity_Window_%04d.png', index));
                export_comparison(connectivity.plv_3d_matrix(:, :, index), ...
                    connectivity.pli_3d_matrix(:, :, index), ...
                    connectivity.wpli_3d_matrix(:, :, index), ...
                    connectivity.channel_names, sprintf( ...
                    'Window %d | center %.4f s', index, ...
                    connectivity.absolute_time_axis(index)), outputFile, ...
                    options.resolution_dpi);
                exportedFiles{end + 1} = outputFile; %#ok<AGROW>
            end
        end
    end

    if options.export_psd && isfield(bundle, 'psd') && ~isempty(bundle.psd)
        source = validatestring(options.psd_source, {'Raw','Filtered'});
        outputFile = fullfile(outputFolder, sprintf('%s_PSD_%s.png', ...
            sanitize(bundle.psd.channel_name), source));
        export_psd(bundle.psd, source, options.psd_in_db, ...
            options.psd_frequency_range_hz, outputFile, options.resolution_dpi);
        exportedFiles{end + 1} = outputFile; %#ok<AGROW>
    end

    if options.export_topoplot && isfield(bundle, 'connectivity') && ...
            ~isempty(bundle.connectivity) && isfield(bundle, 'electrode_mapping') && ...
            bundle.electrode_mapping.num_mapped >= 2
        thresholds = options.topoplot_threshold;
        if options.export_threshold_sweep
            thresholds = 0:0.1:0.9;
        end
        for threshold = thresholds
            outputFile = export_topoplot(bundle, options, threshold, ...
                outputFolder, options.resolution_dpi);
            exportedFiles{end + 1} = outputFile; %#ok<AGROW>
        end
    end

    if options.export_psd_topography && isfield(bundle, 'psd_topography') && ...
            ~isempty(bundle.psd_topography) && ...
            isfield(bundle, 'psd_topography_mapping') && ...
            ~isempty(bundle.psd_topography_mapping) && ...
            bundle.psd_topography_mapping.num_mapped >= 3
        outputFile = export_psd_topography(bundle.psd_topography, ...
            bundle.psd_topography_mapping, options, outputFolder, ...
            options.resolution_dpi);
        exportedFiles{end + 1} = outputFile; %#ok<AGROW>
        if options.export_all_psd_topography_epochs
            outputFile = export_psd_topography_epochs(bundle.psd_topography, ...
                bundle.psd_topography_mapping, options, outputFolder, ...
                options.resolution_dpi);
            exportedFiles{end + 1} = outputFile; %#ok<AGROW>
        end
    end
end

function options = normalize_options(options)
    defaults = struct('export_connectivity', true, 'export_all_windows', false, ...
        'export_psd', true, 'psd_source', 'Raw', 'psd_in_db', false, ...
        'psd_frequency_range_hz', [], ...
        'export_topoplot', true, 'export_threshold_sweep', false, ...
        'topoplot_metric', 'PLV', 'topoplot_scope', 'Whole', ...
        'topoplot_window_index', 1, 'topoplot_threshold', 0.5, ...
        'export_psd_topography', true, ...
        'export_all_psd_topography_epochs', false, ...
        'psd_topography_source', 'Raw', 'psd_topography_in_db', false, ...
        'psd_topography_scope', 'Whole', 'psd_topography_epoch_index', 1, ...
        'resolution_dpi', 300);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index}) || isempty(options.(names{index}))
            options.(names{index}) = defaults.(names{index});
        end
    end
end

function outputFile = export_psd_topography(result, mapping, options, outputFolder, dpi)
    [epochs, whole, mapLimits, colorbarLabel, scaleName] = ...
        psd_topography_display_data(result, mapping, options);
    scope = validatestring(options.psd_topography_scope, {'Whole','Epoch'});
    if strcmp(scope, 'Whole')
        values = whole;
        scopeText = sprintf('Whole %.3f-%.3f s', ...
            result.start_time_sec, result.end_time_sec);
        scopeFile = sprintf('Whole_%.3f-%.3fs', ...
            result.start_time_sec, result.end_time_sec);
    else
        index = min(result.num_epochs, max(1, ...
            round(options.psd_topography_epoch_index)));
        values = epochs(index, :);
        scopeText = sprintf('Epoch %d | %.3f-%.3f s', index, ...
            result.epoch_start_time_sec(index), result.epoch_end_time_sec(index));
        scopeFile = sprintf('Epoch_%04d_%.3f-%.3fs', index, ...
            result.epoch_start_time_sec(index), result.epoch_end_time_sec(index));
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 850 760]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig);
    plot_psd_topography(values, mapping, ax, struct( ...
        'map_limits', mapLimits, 'colorbar_label', colorbarLabel));
    title(ax, sprintf('%s PSD topoplot | %.3g-%.3g Hz | %s | %s', ...
        options.psd_topography_source, result.band_hz(1), result.band_hz(2), ...
        scaleName, scopeText), 'Interpreter', 'none');
    bandName = band_filename(result.band_hz);
    outputFile = fullfile(outputFolder, sprintf( ...
        'PSD_Topoplot_%s_%s_%s_%s.png', ...
        options.psd_topography_source, bandName, scaleName, sanitize(scopeFile)));
    exportgraphics(fig, outputFile, 'Resolution', dpi);
end

function outputFile = export_psd_topography_epochs(result, mapping, options, outputFolder, dpi)
    [epochs, ~, mapLimits, colorbarLabel, scaleName] = ...
        psd_topography_display_data(result, mapping, options);
    count = result.num_epochs;
    columns = min(4, max(1, ceil(sqrt(count))));
    rows = ceil(count / columns);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [50 50 max(900, columns * 330) max(700, rows * 290)]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    layout = tiledlayout(fig, rows, columns, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    lastAxes = [];
    for index = 1:count
        lastAxes = nexttile(layout);
        plot_psd_topography(epochs(index, :), mapping, lastAxes, struct( ...
            'map_limits', mapLimits, 'show_labels', false, 'show_colorbar', false));
        title(lastAxes, sprintf('%d | %.3f-%.3f s', index, ...
            result.epoch_start_time_sec(index), result.epoch_end_time_sec(index)), ...
            'Interpreter', 'none', 'FontSize', 9);
    end
    colorbarHandle = colorbar(lastAxes, 'eastoutside');
    colorbarHandle.Label.String = colorbarLabel;
    title(layout, sprintf('%s PSD topography | %.3g-%.3g Hz | %s', ...
        options.psd_topography_source, result.band_hz(1), result.band_hz(2), ...
        scaleName), 'Interpreter', 'none');
    outputFile = fullfile(outputFolder, sprintf( ...
        'PSD_Topoplot_%s_%s_%s_AllEpochs.png', ...
        options.psd_topography_source, band_filename(result.band_hz), scaleName));
    exportgraphics(fig, outputFile, 'Resolution', dpi);
end

function [epochs, whole, mapLimits, colorbarLabel, scaleName] = ...
        psd_topography_display_data(result, mapping, options)
    source = validatestring(options.psd_topography_source, {'Raw','Filtered'});
    if strcmp(source, 'Filtered')
        epochs = result.filtered_band_power;
        whole = result.filtered_whole_band_power;
    else
        epochs = result.raw_band_power;
        whole = result.raw_whole_band_power;
    end
    if options.psd_topography_in_db
        epochs = 10 * log10(max(epochs, realmin('double')));
        whole = 10 * log10(max(whole, realmin('double')));
        colorbarLabel = 'Integrated band power (dB)';
        scaleName = 'dB';
    else
        colorbarLabel = 'Integrated band power (signal-unit^2)';
        scaleName = 'Linear';
    end
    displayed = [whole; epochs];
    displayed = displayed(:, mapping.analysis_indices);
    mapLimits = [min(displayed(:)), max(displayed(:))];
end

function value = band_filename(band)
    value = sprintf('%g-%gHz', band(1), band(2));
    value = strrep(value, '.', 'p');
end

function export_comparison(plv, pli, wpli, channelNames, titleText, outputFile, dpi)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 100 2100 720]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    matrices = {plv, pli, wpli};
    names = {'PLV','PLI','wPLI'};
    for index = 1:3
        matrix = matrices{index};
        matrix(triu(true(size(matrix)), 1)) = NaN;
        ax = nexttile(layout);
        imageHandle = imagesc(ax, matrix, [0 1]);
        imageHandle.AlphaData = ~isnan(matrix);
        ax.Color = [1 1 1];
        axis(ax, 'square');
        ax.YDir = 'normal';
        ax.XTick = 1:numel(channelNames);
        ax.YTick = 1:numel(channelNames);
        ax.XTickLabel = channelNames;
        ax.YTickLabel = channelNames;
        ax.XTickLabelRotation = 45;
        ax.FontSize = 7;
        title(ax, names{index});
        colorbar(ax);
    end
    colormap(fig, jet(256));
    title(layout, titleText, 'Interpreter', 'none');
    exportgraphics(fig, outputFile, 'Resolution', dpi);
end

function export_psd(result, source, inDb, frequencyRange, outputFile, dpi)
    if strcmp(source, 'Filtered')
        values = result.filtered_psd;
    else
        values = result.raw_psd;
    end
    if inDb
        values = 10 * log10(values + eps);
        unit = 'PSD (dB/Hz)';
    else
        unit = 'PSD (\muV^2/Hz)';
    end
    if isempty(frequencyRange)
        frequencyRange = [result.frequency(1), result.frequency(end)];
    end
    validateattributes(frequencyRange, {'numeric'}, ...
        {'real','finite','vector','numel',2,'nonnegative'}, mfilename, ...
        'psd_frequency_range_hz');
    frequencyRange = double(frequencyRange(:).');
    if frequencyRange(1) >= frequencyRange(2) || ...
            frequencyRange(2) > result.Fs / 2
        error('export_phase_connectivity_results:InvalidPSDFrequencyRange', ...
            'PSD display range must satisfy 0 <= low < high <= Nyquist.');
    end
    visible = result.frequency >= frequencyRange(1) & ...
        result.frequency <= frequencyRange(2);
    if nnz(visible) < 2
        error('export_phase_connectivity_results:EmptyPSDFrequencyRange', ...
            'PSD display range contains too few frequency bins.');
    end
    fig = figure('Visible', 'off', 'Color', 'w');
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig);
    plot(ax, result.frequency(visible), values(visible), 'LineWidth', 1.3);
    grid(ax, 'on');
    xlim(ax, frequencyRange);
    ax.YLimMode = 'auto';
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, unit);
    title(ax, sprintf('%s %s PSD | %.3f-%.3f s', result.channel_name, ...
        source, result.start_time_sec, result.end_time_sec), 'Interpreter', 'none');
    exportgraphics(fig, outputFile, 'Resolution', dpi);
end

function outputFile = export_topoplot(bundle, options, threshold, outputFolder, dpi)
    [matrix, scopeText] = select_matrix(bundle.connectivity, ...
        options.topoplot_metric, options.topoplot_scope, ...
        options.topoplot_window_index);
    mapping = bundle.electrode_mapping;
    mappedMatrix = matrix(mapping.analysis_indices, mapping.analysis_indices);
    mask = tril(true(size(mappedMatrix)), -1) & mappedMatrix >= threshold & mappedMatrix > 0;
    [first, second] = find(mask);
    displayStruct = struct('chanPairs', [first second], ...
        'connectStrength', mappedMatrix(mask), 'connectStrengthLimits', [0 1]);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 850 760]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig);
    stats = topoplot_connect(displayStruct, mapping, ax);
    title(ax, sprintf('%s %s | threshold %.1f | %d connections', ...
        options.topoplot_metric, scopeText, threshold, stats.numConnections), ...
        'Interpreter', 'none');
    scopeName = regexprep(scopeText, '[^A-Za-z0-9]+', '_');
    outputFile = fullfile(outputFolder, sprintf('%s_Topoplot_%s_Threshold_%.1f.png', ...
        options.topoplot_metric, scopeName, threshold));
    exportgraphics(fig, outputFile, 'Resolution', dpi);
end

function [matrix, textValue] = select_matrix(connectivity, metric, scope, windowIndex)
    metric = validatestring(metric, {'PLV','PLI','wPLI'});
    scope = validatestring(scope, {'Whole','Sliding'});
    prefix = lower(metric);
    if strcmpi(metric, 'wPLI')
        prefix = 'wpli';
    end
    if strcmp(scope, 'Whole')
        matrix = connectivity.([prefix '_matrix']);
        textValue = sprintf('Whole_%.3f-%.3fs', connectivity.start_time_sec, ...
            connectivity.end_time_sec);
    else
        count = numel(connectivity.time_axis);
        windowIndex = min(count, max(1, round(windowIndex)));
        values = connectivity.([prefix '_3d_matrix']);
        matrix = values(:, :, windowIndex);
        textValue = sprintf('Window_%04d_%.4fs', windowIndex, ...
            connectivity.absolute_time_axis(windowIndex));
    end
end

function value = sanitize(value)
    value = regexprep(char(string(value)), '[^A-Za-z0-9_-]+', '_');
end
