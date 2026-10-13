function mapping = resolve_electrode_locations(channelNames, source)
%RESOLVE_ELECTRODE_LOCATIONS Map channels to 10-20 or imported locations.
%   MAPPING contains labels, x, y, analysis_indices, unmatched_channels,
%   and source. SOURCE can be empty, a chanlocs structure, or a MAT/SET
%   filename containing chanlocs or EEG.chanlocs.

    channelNames = cellstr(string(channelNames));
    if nargin < 2
        source = [];
    end
    [imported, sourceDescription] = read_source(source);
    builtin = builtin_locations();

    labels = {};
    x = [];
    y = [];
    indices = [];
    unmatched = {};
    for index = 1:numel(channelNames)
        canonical = canonical_label(channelNames{index});
        [found, px, py] = find_location(imported, canonical);
        if ~found
            [found, px, py] = find_location(builtin, canonical);
        end
        if found
            labels{end + 1, 1} = channelNames{index}; %#ok<AGROW>
            x(end + 1, 1) = px; %#ok<AGROW>
            y(end + 1, 1) = py; %#ok<AGROW>
            indices(end + 1, 1) = index; %#ok<AGROW>
        else
            unmatched{end + 1, 1} = channelNames{index}; %#ok<AGROW>
        end
    end

    mapping = struct('labels', {labels}, 'x', x, 'y', y, ...
        'analysis_indices', indices, 'unmatched_channels', {unmatched}, ...
        'source', sourceDescription, 'num_mapped', numel(indices));
end

function [locations, description] = read_source(source)
    locations = struct('labels', {}, 'x', {}, 'y', {});
    description = 'Built-in 10-20 montage';
    if isempty(source)
        return;
    end
    if isstruct(source)
        chanlocs = source;
        description = 'Imported chanlocs structure';
    else
        fileName = char(string(source));
        if ~isfile(fileName)
            error('resolve_electrode_locations:FileNotFound', ...
                'Montage file was not found: %s', fileName);
        end
        [~, ~, extension] = fileparts(fileName);
        if ~ismember(lower(extension), {'.mat','.set'})
            error('resolve_electrode_locations:UnsupportedFileType', ...
                'Montage source must be a MAT or SET file.');
        end
        variableInfo = whos('-file', fileName);
        variableNames = {variableInfo.name};
        if any(strcmp(variableNames, 'chanlocs'))
            fileObject = matfile(fileName);
            chanlocs = fileObject.chanlocs;
        elseif any(strcmp(variableNames, 'EEG'))
            loaded = load(fileName, 'EEG');
            if ~isfield(loaded.EEG, 'chanlocs')
                error('resolve_electrode_locations:MissingChanlocs', ...
                    'EEG structure does not contain chanlocs.');
            end
            chanlocs = loaded.EEG.chanlocs;
        else
            error('resolve_electrode_locations:MissingChanlocs', ...
                'Montage file must contain chanlocs or EEG.chanlocs.');
        end
        description = fileName;
    end
    locations = convert_chanlocs(chanlocs);
end

function locations = convert_chanlocs(chanlocs)
    if ~isstruct(chanlocs) || isempty(chanlocs) || ~isfield(chanlocs, 'labels')
        error('resolve_electrode_locations:InvalidChanlocs', ...
            'chanlocs must be a nonempty structure with labels.');
    end
    locations = repmat(struct('labels', '', 'x', NaN, 'y', NaN), numel(chanlocs), 1);
    for index = 1:numel(chanlocs)
        locations(index).labels = char(string(chanlocs(index).labels));
        if has_finite_fields(chanlocs(index), {'theta','radius'})
            theta = double(chanlocs(index).theta) * pi / 180;
            radius = double(chanlocs(index).radius);
            locations(index).x = sin(theta) * radius;
            locations(index).y = cos(theta) * radius;
        elseif has_finite_fields(chanlocs(index), {'X','Y'})
            locations(index).x = -double(chanlocs(index).Y);
            locations(index).y = double(chanlocs(index).X);
        else
            continue;
        end
    end
    valid = isfinite([locations.x]) & isfinite([locations.y]);
    locations = locations(valid);
    if isempty(locations)
        error('resolve_electrode_locations:NoCoordinates', ...
            'No finite theta/radius or X/Y coordinates were found.');
    end
    radial = hypot([locations.x], [locations.y]);
    scale = max(radial);
    if scale > 0
        for index = 1:numel(locations)
            locations(index).x = locations(index).x / scale;
            locations(index).y = locations(index).y / scale;
        end
    end
end

function tf = has_finite_fields(item, names)
    tf = all(isfield(item, names));
    if ~tf
        return;
    end
    for index = 1:numel(names)
        value = item.(names{index});
        tf = tf && isnumeric(value) && isscalar(value) && isfinite(value);
    end
end

function [found, x, y] = find_location(locations, canonical)
    found = false;
    x = NaN;
    y = NaN;
    for index = 1:numel(locations)
        if strcmp(canonical_label(locations(index).labels), canonical)
            found = true;
            x = locations(index).x;
            y = locations(index).y;
            return;
        end
    end
end

function locations = builtin_locations()
    labels = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','Cz', ...
        'C4','T4','T5','P3','Pz','P4','T6','O1','O2'};
    xy = [-.30 .95; .30 .95; -.85 .55; -.45 .52; 0 .58; .45 .52; .85 .55; ...
        -1 0; -.50 0; 0 0; .50 0; 1 0; -.85 -.52; -.45 -.52; 0 -.58; ...
        .45 -.52; .85 -.52; -.30 -.95; .30 -.95];
    locations = repmat(struct('labels', '', 'x', 0, 'y', 0), numel(labels), 1);
    for index = 1:numel(labels)
        locations(index).labels = labels{index};
        locations(index).x = xy(index, 1);
        locations(index).y = xy(index, 2);
    end
end

function value = canonical_label(value)
    value = upper(regexprep(char(string(value)), '[^A-Za-z0-9]', ''));
    value = regexprep(value, 'AVG$', '');
end
