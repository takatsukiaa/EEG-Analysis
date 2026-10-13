%% PhaseConnectivity_runner
% Calculate PLV, PLI, and wPLI from the same filtered EEG samples, Hilbert
% analytic signals, and sliding windows for direct comparison.

%% Validated-backend compatibility wrapper
% Set PHASE_CONNECTIVITY_USE_LEGACY_SCRIPT=true before running to execute
% the preserved standalone implementation below.
if ~exist('PHASE_CONNECTIVITY_USE_LEGACY_SCRIPT', 'var') || ...
        ~PHASE_CONNECTIVITY_USE_LEGACY_SCRIPT
    clc; close all;
    script_dir = fileparts(mfilename('fullpath'));
    input_file = fullfile(script_dir, 'Raw_Data', 'No1.edf');
    recording = load_phase_connectivity_recording(input_file);
    options = struct('start_time_sec', 144.5, 'end_time_sec', 146, ...
        'channel_names', {recording.channel_names}, ...
        'window_length_sec', 0.5, 'window_step_sec', 0.1, ...
        'filter_settings', struct('bandpass_enabled', true, ...
        'band_hz', [25 50], 'bandpass_order', 20, 'notch_hz', 0));
    connectivity_results = run_phase_connectivity_analysis(recording, options);
    plv_matrix = connectivity_results.plv_matrix;
    pli_matrix = connectivity_results.pli_matrix;
    wpli_matrix = connectivity_results.wpli_matrix;
    plv_3d_matrix = connectivity_results.plv_3d_matrix;
    pli_3d_matrix = connectivity_results.pli_3d_matrix;
    wpli_3d_matrix = connectivity_results.wpli_3d_matrix;
    t_axis = connectivity_results.time_axis;
    channel_names = connectivity_results.channel_names;
    fprintf('Aligned phase-connectivity run complete: %d channels, %d windows.\n', ...
        numel(channel_names), numel(t_axis));
    return;
end

clc;
clearvars;
close all;

script_full_path = mfilename('fullpath');
if isempty(script_full_path)
    script_dir = pwd;
else
    script_dir = fileparts(script_full_path);
end

%% Shared configuration
input_file = fullfile('Raw_Data', 'No1.edf');
start_time_sec = 144.5;
end_time_sec = 146;

frequency_band_hz = [25 50];
filter_order = 20;

window_length_sec = 0.5;
window_step_sec = 0.1;

% Comparison figures are enabled by default. Set the environment variable
% PHASE_CONNECTIVITY_PLOT_RESULTS=false to suppress file export.
plot_results = ~strcmpi(getenv('PHASE_CONNECTIVITY_PLOT_RESULTS'), 'false');
plot_resolution_dpi = 300;

%% Resolve and load the EEG file
if ~isfile(input_file)
    input_file = fullfile(script_dir, input_file);
end
if ~isfile(input_file)
    error('PhaseConnectivity_runner:FileNotFound', ...
        'EEG file not found: %s', input_file);
end

[eeg_data, Fs] = load_eeg_file(input_file);
eeg_data = order_eeg_channels(eeg_data);
[num_recording_samples, channel_names] = validate_eeg_structure(eeg_data);

fprintf('Loaded %s\n', input_file);
fprintf('Sampling rate: %.6g Hz; EEG channels: %d; samples: %d\n', ...
    Fs, numel(channel_names), num_recording_samples);

%% Validate shared parameters and extract [start, end) samples
validateattributes(start_time_sec, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'start_time_sec');
validateattributes(end_time_sec, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'end_time_sec');
if end_time_sec <= start_time_sec
    error('PhaseConnectivity_runner:InvalidSegment', ...
        'end_time_sec (%.6g) must be greater than start_time_sec (%.6g).', ...
        end_time_sec, start_time_sec);
end

recording_duration_sec = num_recording_samples / Fs;
if end_time_sec > recording_duration_sec
    error('PhaseConnectivity_runner:SegmentExceedsRecording', ...
        ['The requested end time (%.6g s) exceeds the recording duration ' ...
         '(%.6g s).'], end_time_sec, recording_duration_sec);
end

validateattributes(frequency_band_hz, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', 2, 'positive'}, ...
    mfilename, 'frequency_band_hz');
frequency_band_hz = double(frequency_band_hz(:).');
if frequency_band_hz(1) >= frequency_band_hz(2)
    error('PhaseConnectivity_runner:InvalidFrequencyBand', ...
        'frequency_band_hz must be [low high] with low < high.');
end
if frequency_band_hz(2) >= Fs / 2
    error('PhaseConnectivity_runner:FrequencyExceedsNyquist', ...
        ['The upper band edge (%.6g Hz) must be below the Nyquist ' ...
         'frequency (%.6g Hz).'], frequency_band_hz(2), Fs / 2);
end

validateattributes(filter_order, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'filter_order');
validateattributes(window_length_sec, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'window_length_sec');
validateattributes(window_step_sec, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'window_step_sec');

start_sample = floor(start_time_sec * Fs) + 1;
end_sample = min(floor(end_time_sec * Fs), num_recording_samples);
segment_samples = end_sample - start_sample + 1;
if segment_samples < 2
    error('PhaseConnectivity_runner:SegmentTooShort', ...
        'The selected interval contains fewer than two samples.');
end
if segment_samples <= 3 * filter_order
    error('PhaseConnectivity_runner:SegmentTooShortForFilter', ...
        ['The selected segment has %d samples. Zero-phase filtering with ' ...
         'order %d requires more than approximately %d samples.'], ...
        segment_samples, filter_order, 3 * filter_order);
end

eeg_segment = struct();
for channel_index = 1:numel(channel_names)
    channel_name = channel_names{channel_index};
    signal = eeg_data.(channel_name);
    eeg_segment.(channel_name) = double(signal(start_sample:end_sample));
end

fprintf('Shared segment: %.6g-%.6g s (%d samples)\n', ...
    start_time_sec, end_time_sec, segment_samples);
fprintf('Shared band-pass filter: %.6g-%.6g Hz, FIR order %d\n', ...
    frequency_band_hz(1), frequency_band_hz(2), filter_order);

%% Filter once and calculate all metrics together
filtered_eeg = band_pass(eeg_segment, frequency_band_hz(1), ...
    frequency_band_hz(2), Fs, filter_order);

connectivity_results = get_phase_connectivity( ...
    filtered_eeg, Fs, window_length_sec, window_step_sec);

plv_matrix = connectivity_results.plv_matrix;
pli_matrix = connectivity_results.pli_matrix;
wpli_matrix = connectivity_results.wpli_matrix;
plv_3d_matrix = connectivity_results.plv_3d_matrix;
pli_3d_matrix = connectivity_results.pli_3d_matrix;
wpli_3d_matrix = connectivity_results.wpli_3d_matrix;
t_axis = connectivity_results.time_axis;
channel_names = connectivity_results.channel_names;

fprintf('Whole matrices: %d x %d for PLV, PLI, and wPLI\n', ...
    size(plv_matrix, 1), size(plv_matrix, 2));
fprintf('Aligned sliding matrices: %d windows (%.6g s window, %.6g s step)\n', ...
    numel(t_axis), window_length_sec, window_step_sec);

%% Optional three-panel comparison figures
if plot_results
    [~, input_stem] = fileparts(input_file);
    analysis_tag = sprintf('%s_%g-%gs_%g-%gHz', input_stem, ...
        start_time_sec, end_time_sec, frequency_band_hz(1), frequency_band_hz(2));
    analysis_tag = regexprep(analysis_tag, '[^A-Za-z0-9_-]', '_');
    output_folder = fullfile(script_dir, 'PhaseConnectivity_Results', analysis_tag);
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    whole_title = sprintf( ...
        'Whole segment %.3f-%.3f s | %.3f-%.3f Hz', ...
        start_time_sec, end_time_sec, frequency_band_hz(1), frequency_band_hz(2));
    export_comparison_figure(plv_matrix, pli_matrix, wpli_matrix, ...
        channel_names, whole_title, ...
        fullfile(output_folder, 'PhaseConnectivity_Whole.png'), ...
        plot_resolution_dpi);

    for window_index = 1:numel(t_axis)
        absolute_time_sec = start_time_sec + t_axis(window_index);
        window_title = sprintf( ...
            'Window %d | center %.4f s | %.3f-%.3f Hz', ...
            window_index, absolute_time_sec, ...
            frequency_band_hz(1), frequency_band_hz(2));
        export_comparison_figure( ...
            plv_3d_matrix(:, :, window_index), ...
            pli_3d_matrix(:, :, window_index), ...
            wpli_3d_matrix(:, :, window_index), ...
            channel_names, window_title, ...
            fullfile(output_folder, sprintf( ...
                'PhaseConnectivity_Window_%03d.png', window_index)), ...
            plot_resolution_dpi);
    end
    fprintf('Comparison figures saved to: %s\n', output_folder);
end

fprintf(['PLV/PLI/wPLI comparison complete. Aligned matrices and t_axis ' ...
    'remain in the MATLAB workspace.\n']);

%% Local file-loading, validation, and plotting helpers
function [eeg_data, Fs] = load_eeg_file(input_file)
    [~, ~, extension] = fileparts(input_file);

    switch lower(extension)
        case '.mat'
            loaded = load(input_file);
            [Fs, sample_rate_field] = read_mat_sample_rate(loaded);

            excluded_fields = {sample_rate_field, 'data', ...
                'C4_segment1', 'C4_segment2', 'O1_segment1', 'O1_segment2', ...
                'O2_segment1', 'O2_segment2', 'P4_segment1', 'P4_segment2'};
            eeg_data = retain_eeg_vector_fields(loaded, excluded_fields);

        case '.edf'
            info = edfinfo(input_file);
            record_duration_sec = seconds(info.DataRecordDuration);
            if ~isscalar(record_duration_sec) || ~isfinite(record_duration_sec) || ...
                    record_duration_sec <= 0
                error('PhaseConnectivity_runner:InvalidEDFRecordDuration', ...
                    'EDF data-record duration must be a positive scalar.');
            end

            sample_rates = double(info.NumSamples(:)) / record_duration_sec;
            edf_table = edfread(input_file);
            variable_names = edf_table.Properties.VariableNames;
            eeg_data = struct();
            retained_rates = [];

            for variable_index = 1:numel(variable_names)
                channel_name = variable_names{variable_index};
                if is_non_eeg_name(channel_name) || variable_index > numel(sample_rates)
                    continue;
                end

                column = edf_table.(channel_name);
                if iscell(column)
                    if isempty(column) || any(~cellfun(@isnumeric, column))
                        continue;
                    end
                    blocks = cellfun(@(x) double(x(:)), column, ...
                        'UniformOutput', false);
                    signal = vertcat(blocks{:});
                elseif isnumeric(column)
                    signal = double(column(:));
                else
                    continue;
                end

                if isempty(signal) || ~isreal(signal) || any(~isfinite(signal))
                    error('PhaseConnectivity_runner:InvalidEDFChannel', ...
                        'EDF channel "%s" is empty, complex, or nonfinite.', channel_name);
                end
                eeg_data.(channel_name) = signal;
                retained_rates(end + 1) = sample_rates(variable_index); %#ok<AGROW>
            end

            if isempty(retained_rates)
                error('PhaseConnectivity_runner:NoEEGChannels', ...
                    'No numeric EEG channels were found in the EDF file.');
            end
            Fs = retained_rates(1);
            rate_tolerance = max(1, abs(Fs)) * 1e-10;
            if any(abs(retained_rates - Fs) > rate_tolerance)
                error('PhaseConnectivity_runner:MixedSamplingRates', ...
                    ['Retained EDF channels have different sampling rates. ' ...
                     'Resample them to a common rate before calculating connectivity.']);
            end

        otherwise
            error('PhaseConnectivity_runner:UnsupportedFileType', ...
                'Unsupported EEG file extension "%s". Use .edf or .mat.', extension);
    end

    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'Fs');
end

function [Fs, sample_rate_field] = read_mat_sample_rate(loaded)
    if isfield(loaded, 'Fs')
        sample_rate_field = 'Fs';
    elseif isfield(loaded, 'fs')
        sample_rate_field = 'fs';
    else
        error('PhaseConnectivity_runner:MissingSamplingRate', ...
            'MAT file must contain a scalar sampling rate named Fs or fs.');
    end
    Fs = loaded.(sample_rate_field);
    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, sample_rate_field);
    Fs = double(Fs);
end

function eeg_data = retain_eeg_vector_fields(loaded, excluded_fields)
    eeg_data = struct();
    all_fields = fieldnames(loaded);
    for field_index = 1:numel(all_fields)
        field_name = all_fields{field_index};
        if any(strcmpi(field_name, excluded_fields)) || is_non_eeg_name(field_name)
            continue;
        end

        value = loaded.(field_name);
        if isnumeric(value) && isreal(value) && isvector(value) && ~isempty(value)
            if any(~isfinite(value(:)))
                error('PhaseConnectivity_runner:NonfiniteMATChannel', ...
                    'MAT channel "%s" contains NaN or Inf values.', field_name);
            end
            eeg_data.(field_name) = double(value(:));
        end
    end
    if isempty(fieldnames(eeg_data))
        error('PhaseConnectivity_runner:NoEEGChannels', ...
            'No numeric EEG channel vectors were found in the MAT file.');
    end
end

function tf = is_non_eeg_name(channel_name)
    tf = contains(channel_name, 'emg', 'IgnoreCase', true) || ...
         contains(channel_name, 'ekg', 'IgnoreCase', true) || ...
         contains(channel_name, 'ecg', 'IgnoreCase', true) || ...
         contains(channel_name, 'annotation', 'IgnoreCase', true);
end

function eeg_data = order_eeg_channels(eeg_data)
    preferred_order = {'Fp1_AVG', 'Fp2_AVG', 'F7_AVG', 'F3_AVG', 'FZ_AVG', ...
        'F4_AVG', 'F8_AVG', 'T3_AVG', 'C3_AVG', 'CZ_AVG', 'C4_AVG', ...
        'T4_AVG', 'T5_AVG', 'P3_AVG', 'PZ_AVG', 'P4_AVG', 'T6_AVG', ...
        'O1_AVG', 'O2_AVG'};
    current_names = fieldnames(eeg_data);
    used = false(size(current_names));
    ordered_names = cell(0, 1);

    for preferred_index = 1:numel(preferred_order)
        match = find(strcmpi(current_names, preferred_order{preferred_index}) & ...
            ~used, 1);
        if ~isempty(match)
            ordered_names{end + 1, 1} = current_names{match}; %#ok<AGROW>
            used(match) = true;
        end
    end
    ordered_names = [ordered_names; current_names(~used)];
    eeg_data = orderfields(eeg_data, ordered_names);
end

function [num_samples, channel_names] = validate_eeg_structure(eeg_data)
    channel_names = fieldnames(eeg_data);
    if numel(channel_names) < 2
        error('PhaseConnectivity_runner:NotEnoughChannels', ...
            'At least two EEG channels are required; found %d.', numel(channel_names));
    end

    num_samples = numel(eeg_data.(channel_names{1}));
    if num_samples < 2
        error('PhaseConnectivity_runner:NotEnoughSamples', ...
            'EEG channels must contain at least two samples.');
    end
    for channel_index = 1:numel(channel_names)
        signal = eeg_data.(channel_names{channel_index});
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('PhaseConnectivity_runner:InvalidChannel', ...
                'Channel "%s" must be a finite, real numeric vector.', ...
                channel_names{channel_index});
        end
        if numel(signal) ~= num_samples
            error('PhaseConnectivity_runner:UnequalChannelLengths', ...
                'Channel "%s" has %d samples; expected %d.', ...
                channel_names{channel_index}, numel(signal), num_samples);
        end
    end
end

function export_comparison_figure(plv, pli, wpli, channel_names, ...
        title_text, output_file, resolution_dpi)
    matrices = {plv, pli, wpli};
    metric_names = {'PLV', 'PLI', 'wPLI'};
    upper_mask = triu(true(size(plv)), 1);

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [50 100 2100 720]);
    close_figure = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for metric_index = 1:numel(matrices)
        displayed_matrix = matrices{metric_index};
        displayed_matrix(upper_mask) = NaN;

        ax = nexttile(layout, metric_index);
        image_handle = imagesc(ax, displayed_matrix, [0 1]);
        image_handle.AlphaData = ~isnan(displayed_matrix);
        ax.Color = [1 1 1];
        axis(ax, 'square');
        ax.YDir = 'normal';
        ax.XTick = 1:numel(channel_names);
        ax.YTick = 1:numel(channel_names);
        ax.XTickLabel = channel_names;
        ax.YTickLabel = channel_names;
        ax.XTickLabelRotation = 45;
        ax.FontSize = 7;
        title(ax, metric_names{metric_index});
        colorbar(ax);
    end

    colormap(fig, jet(256));
    title(layout, title_text, 'Interpreter', 'none');
    exportgraphics(fig, output_file, 'Resolution', resolution_dpi);
end
