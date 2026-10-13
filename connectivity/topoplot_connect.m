function stats = topoplot_connect(displayStruct, locations, parentAxes)
%TOPOPLOT_CONNECT Draw a thresholded connectivity network on a scalp.
%   STATS = TOPOPLOT_CONNECT(DISPLAYSTRUCT, LOCATIONS, PARENTAXES) renders
%   channel pairs on the specified axes without requiring EEGLAB. LOCATIONS
%   can be a mapping returned by RESOLVE_ELECTRODE_LOCATIONS or an EEGLAB
%   chanlocs structure. The two-argument form draws into the current axes.

    if nargin < 3 || isempty(parentAxes)
        parentAxes = gca;
    end
    [labels, x, y] = normalize_locations(locations);

    if ~isstruct(displayStruct) || ~isfield(displayStruct, 'chanPairs')
        error('topoplot_connect:InvalidDisplayStruct', ...
            'displayStruct must contain chanPairs.');
    end
    pairs = displayStruct.chanPairs;
    if isempty(pairs)
        pairs = zeros(0, 2);
    end
    if ~isnumeric(pairs) || size(pairs, 2) ~= 2 || ...
            any(pairs(:) < 1) || any(pairs(:) > numel(x)) || ...
            any(fix(pairs(:)) ~= pairs(:))
        error('topoplot_connect:InvalidPairs', ...
            'chanPairs must contain valid one-based channel index pairs.');
    end
    if isfield(displayStruct, 'connectStrength')
        strengths = double(displayStruct.connectStrength(:));
    else
        strengths = 0.5 * ones(size(pairs, 1), 1);
    end
    if numel(strengths) ~= size(pairs, 1) || any(~isfinite(strengths))
        error('topoplot_connect:InvalidStrengths', ...
            'connectStrength must contain one finite value per channel pair.');
    end
    if isfield(displayStruct, 'connectStrengthLimits')
        limits = double(displayStruct.connectStrengthLimits(:).');
    else
        limits = [0 1];
    end
    if numel(limits) ~= 2 || limits(1) >= limits(2)
        error('topoplot_connect:InvalidLimits', ...
            'connectStrengthLimits must be an increasing two-element vector.');
    end

    cla(parentAxes);
    hold(parentAxes, 'on');
    parentAxes.Color = [1 1 1];
    parentAxes.Visible = 'off';
    axis(parentAxes, 'equal');
    xlim(parentAxes, [-1.28 1.28]);
    ylim(parentAxes, [-1.18 1.28]);

    theta = linspace(0, 2*pi, 361);
    plot(parentAxes, cos(theta), sin(theta), 'k-', 'LineWidth', 1.8);
    plot(parentAxes, [-.18 0 .18], [.98 1.20 .98], 'k-', 'LineWidth', 1.8);
    plot(parentAxes, [-1.00 -1.10 -1.12 -1.04 -1.00], ...
        [.35 .28 0 -.28 -.35], 'k-', 'LineWidth', 1.5);
    plot(parentAxes, [1.00 1.10 1.12 1.04 1.00], ...
        [.35 .28 0 -.28 -.35], 'k-', 'LineWidth', 1.5);

    colorMap = jet(256);
    for index = 1:size(pairs, 1)
        normalized = (strengths(index) - limits(1)) / diff(limits);
        colorIndex = 1 + round(min(1, max(0, normalized)) * (size(colorMap, 1) - 1));
        plot(parentAxes, x(pairs(index, :)), y(pairs(index, :)), '-', ...
            'LineWidth', 1.5, 'Color', colorMap(colorIndex, :));
    end

    scatter(parentAxes, x, y, 36, 'k', 'filled', ...
        'MarkerFaceColor', [0.08 0.08 0.08]);
    for index = 1:numel(labels)
        text(parentAxes, x(index), y(index) - 0.055, labels{index}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Interpreter', 'none');
    end
    colormap(parentAxes, colorMap);
    clim(parentAxes, limits);
    colorbar(parentAxes);
    hold(parentAxes, 'off');

    if isempty(strengths)
        minimum = NaN;
        maximum = NaN;
    else
        minimum = min(strengths);
        maximum = max(strengths);
    end
    stats = struct('numConnections', numel(strengths), ...
        'minimumStrength', minimum, 'maximumStrength', maximum);
end

function [labels, x, y] = normalize_locations(locations)
    if isstruct(locations) && isscalar(locations) && ...
            all(isfield(locations, {'labels','x','y'})) && ...
            iscell(locations.labels)
        labels = locations.labels;
        x = double(locations.x(:));
        y = double(locations.y(:));
        return;
    end
    if ~isstruct(locations) || isempty(locations) || ~isfield(locations, 'labels')
        error('topoplot_connect:InvalidLocations', ...
            'locations must be a resolved mapping or EEGLAB chanlocs structure.');
    end
    labels = cellstr(string({locations.labels}));
    x = NaN(numel(locations), 1);
    y = NaN(numel(locations), 1);
    for index = 1:numel(locations)
        if all(isfield(locations(index), {'theta','radius'})) && ...
                isfinite_numeric(locations(index).theta) && ...
                isfinite_numeric(locations(index).radius)
            angle = double(locations(index).theta) * pi / 180;
            x(index) = sin(angle) * double(locations(index).radius);
            y(index) = cos(angle) * double(locations(index).radius);
        elseif all(isfield(locations(index), {'X','Y'})) && ...
                isfinite_numeric(locations(index).X) && ...
                isfinite_numeric(locations(index).Y)
            x(index) = -double(locations(index).Y);
            y(index) = double(locations(index).X);
        end
    end
    valid = isfinite(x) & isfinite(y);
    labels = labels(valid);
    x = x(valid);
    y = y(valid);
    scale = max(hypot(x, y));
    if scale > 0
        x = x / scale;
        y = y / scale;
    end
end

function tf = isfinite_numeric(value)
    tf = isnumeric(value) && isscalar(value) && isfinite(value);
end
