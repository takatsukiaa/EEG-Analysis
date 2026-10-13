function stats = plot_psd_topography(values, mapping, parentAxes, options)
%PLOT_PSD_TOPOGRAPHY Render a scalar scalp map without EEGLAB.
%   VALUES are in analysis-channel order. MAPPING is returned by
%   resolve_electrode_locations. PARENTAXES may be a UIAxes or regular axes.

    if nargin < 4 || isempty(options)
        options = struct();
    end
    options = normalize_options(options);
    validateattributes(values, {'numeric'}, {'real','vector','nonempty'}, ...
        mfilename, 'values');
    if ~isstruct(mapping) || ~isscalar(mapping) || ...
            ~all(isfield(mapping, {'x','y','labels','analysis_indices','num_mapped'}))
        error('plot_psd_topography:InvalidMapping', ...
            'mapping must be returned by resolve_electrode_locations.');
    end
    if ~isgraphics(parentAxes, 'axes')
        error('plot_psd_topography:InvalidAxes', ...
            'parentAxes must be a valid axes or UIAxes.');
    end

    values = double(values(:));
    if numel(values) == mapping.num_mapped
        mappedValues = values;
    elseif ~isempty(mapping.analysis_indices) && ...
            numel(values) >= max(mapping.analysis_indices)
        mappedValues = values(mapping.analysis_indices);
    else
        error('plot_psd_topography:ValueCountMismatch', ...
            'values do not match the analysis or mapped channel count.');
    end
    x = double(mapping.x(:));
    y = double(mapping.y(:));
    valid = isfinite(x) & isfinite(y) & isfinite(mappedValues);
    x = x(valid); y = y(valid); mappedValues = mappedValues(valid);
    labels = mapping.labels(valid);
    if numel(mappedValues) < 3
        error('plot_psd_topography:InsufficientChannels', ...
            'At least three mapped channels with finite values are required.');
    end
    if rank([x - mean(x), y - mean(y)]) < 2
        error('plot_psd_topography:CollinearChannels', ...
            'Mapped channels must not all lie on one line.');
    end

    mapLimits = options.map_limits;
    if isempty(mapLimits)
        mapLimits = [min(mappedValues), max(mappedValues)];
    end
    validateattributes(mapLimits, {'numeric'}, ...
        {'real','finite','vector','numel',2}, mfilename, 'map_limits');
    mapLimits = double(mapLimits(:).');
    if mapLimits(1) > mapLimits(2)
        error('plot_psd_topography:InvalidMapLimits', ...
            'map_limits must be increasing.');
    end
    if mapLimits(1) == mapLimits(2)
        padding = max(abs(mapLimits(1)) * 0.01, 1e-12);
        mapLimits = mapLimits + [-padding padding];
    end

    gridAxis = linspace(-1, 1, options.grid_size);
    [gridX, gridY] = meshgrid(gridAxis, gridAxis);
    interpolant = scatteredInterpolant(x, y, mappedValues, 'natural', 'nearest');
    gridValues = interpolant(gridX, gridY);
    gridValues(hypot(gridX, gridY) > 1) = NaN;

    cla(parentAxes, 'reset');
    surface(parentAxes, gridX, gridY, zeros(size(gridValues)), gridValues, ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    view(parentAxes, 2);
    hold(parentAxes, 'on');
    theta = linspace(0, 2*pi, 300);
    plot(parentAxes, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5);
    plot(parentAxes, [-.12 0 .12], [.99 1.14 .99], 'k-', 'LineWidth', 1.5);
    plot(parentAxes, [-1 -1.08 -1.10 -1.08 -1], [.18 .12 0 -.12 -.18], ...
        'k-', 'LineWidth', 1.2);
    plot(parentAxes, [1 1.08 1.10 1.08 1], [.18 .12 0 -.12 -.18], ...
        'k-', 'LineWidth', 1.2);
    scatter(parentAxes, x, y, 24, 'k', 'filled', ...
        'MarkerFaceColor', [0.98 0.98 0.98], 'LineWidth', 0.8);
    if options.show_labels
        text(parentAxes, x, y - 0.045, labels, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none');
    end
    hold(parentAxes, 'off');
    axis(parentAxes, 'equal');
    xlim(parentAxes, [-1.15 1.15]);
    ylim(parentAxes, [-1.12 1.18]);
    parentAxes.Visible = 'off';
    clim(parentAxes, mapLimits);
    colormap(parentAxes, jet(256));
    if options.show_colorbar
        colorbarHandle = colorbar(parentAxes);
        colorbarHandle.Label.String = options.colorbar_label;
    end

    stats = struct('num_mapped', numel(mappedValues), ...
        'minimum_value', min(mappedValues), 'maximum_value', max(mappedValues), ...
        'map_limits', mapLimits);
end

function options = normalize_options(options)
    if ~isstruct(options) || ~isscalar(options)
        error('plot_psd_topography:InvalidOptions', ...
            'options must be a scalar structure.');
    end
    defaults = struct('map_limits', [], 'show_labels', true, ...
        'show_colorbar', true, 'colorbar_label', 'Band power', 'grid_size', 151);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index}) || isempty(options.(names{index}))
            options.(names{index}) = defaults.(names{index});
        end
    end
    validateattributes(options.grid_size, {'numeric'}, ...
        {'real','finite','scalar','integer','>=',25}, mfilename, 'grid_size');
    options.show_labels = logical(options.show_labels);
    options.show_colorbar = logical(options.show_colorbar);
    options.colorbar_label = char(string(options.colorbar_label));
end
