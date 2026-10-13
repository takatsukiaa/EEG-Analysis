%% EEG_Windowed_MTLCCA_ArtifactRemoval_runner
% Windowed multi-time-lag CCA for conservative muscle-artifact attenuation.
% Edit only the configuration block for a different file or interval.
%
% EMG is optional component evidence only; the MCCA decomposition uses EEG.
% During a GTC, fast ictal activity can resemble muscle. Compare every result
% with the untouched recording before clinical interpretation.

clc;
clearvars;
close all;

script_full_path = mfilename('fullpath');
if isempty(script_full_path)
    script_dir = pwd;
else
    script_dir = fileparts(script_full_path);
end

%% Configuration
% Leave empty to select an EDF or MAT file from a dialog. Set an explicit
% path here only when running the script noninteractively.
input_file = getenv('EEG_MTLCCA_INPUT_FILE');

% Half-open interval [start_time_sec, end_time_sec). Set either value to []
% to use the corresponding recording boundary; set both to [] for all data.
start_time_sec = [];
end_time_sec = [];

window_duration_sec = 8;
overlap_fraction = 0.75;
lag_samples = [];                 % [] gives 10 lags spanning 0-40 ms.
num_lag_datasets = 10;
max_lag_duration_sec = 0.040;
covariance_regularization = 1e-6;
mcca_max_sweeps = 200;
mcca_tolerance = 1e-8;
mcca_initializations = 3;
max_candidate_fraction = 0.25;
minimum_component_gain = 0.25;
apply_preservation_guard = true;
minimum_low_band_correlation = 0.98;
maximum_low_band_power_change_db = 0.5;
minimum_spatial_covariance_correlation = 0.95;
robust_clip_mad = 12;

% [] uses automatic gains. A scalar/vector applies globally; a cell array
% supplies one scalar/vector per window. Guardrails remain enabled by default.
component_gains_override = [];

save_results = true;
save_edf_result = true;
plot_results = true;
plot_resolution_dpi = 300;

% Column labels for matrix-based MAT files. This ordering was verified for
% Raw_Data/No2.mat. Edit it when another MAT file uses a different layout.
mat_channel_names = {"Fp1" "Fp2" "Fz" "F7" "F3" "F4" "F8" "T3" "C3" "Cz" "C4" "T4" "T5" "P3" "Pz" "P4" "T6" "O1" "O2", "EMG1 +", "EMG1 -"};

%% Resolve and load the recording
if isempty(input_file)
    initial_folder = fullfile(script_dir, 'Raw_Data');
    if ~isfolder(initial_folder)
        initial_folder = script_dir;
    end
    file_filter = { ...
        '*.edf;*.mat', 'EEG recordings (*.edf, *.mat)'; ...
        '*.edf', 'European Data Format (*.edf)'; ...
        '*.mat', 'MATLAB data files (*.mat)'; ...
        '*.*', 'All files (*.*)'};
    [selected_file, selected_folder] = uigetfile(file_filter, ...
        'Select an EEG recording for windowed multi-time-lag CCA', initial_folder);
    if isequal(selected_file, 0)
        fprintf('File selection cancelled. No data were processed.\n');
        return;
    end
    input_file = fullfile(selected_folder, selected_file);
elseif ~isfile(input_file)
    input_file = fullfile(script_dir, input_file);
end
if ~isfile(input_file)
    error('EEG_BSS_CCA_ArtifactRemoval_runner:FileNotFound', ...
        'Input file not found: %s', input_file);
end

[all_signals, Fs, source_metadata] = load_recording( ...
    input_file, mat_channel_names);
[eeg_data, emg_data, auxiliary_data] = split_recording_channels(all_signals);
eeg_data = order_eeg_channels(eeg_data);
[num_recording_samples, eeg_channel_names] = validate_signal_structure( ...
    eeg_data, 'EEG');

selected_emg_name = '';
selected_emg_signal = zeros(0, 1);
if ~isempty(fieldnames(emg_data))
    try
        [selected_emg_name, selected_emg_signal] = ...
            select_emg_reference(emg_data, '');
        if numel(selected_emg_signal) ~= num_recording_samples
            error('EEG_BSS_CCA_ArtifactRemoval_runner:UnequalEEGAndEMGLengths', ...
                ['EEG channels contain %d samples while the optional EMG ' ...
                 'reference contains %d.'], num_recording_samples, ...
                numel(selected_emg_signal));
        end
    catch exception
        warning('EEG_BSS_CCA_ArtifactRemoval_runner:EMGDiagnosticUnavailable', ...
            ['No unique optional EMG diagnostic reference was selected (%s). ' ...
             'MTLCCA processing will continue because it uses EEG only.'], ...
            exception.message);
        selected_emg_name = '';
        selected_emg_signal = zeros(0, 1);
    end
end

recording_duration_sec = num_recording_samples / Fs;
fprintf('Loaded %s\n', input_file);
fprintf('Sampling rate: %.6g Hz; EEG channels: %d; samples: %d\n', ...
    Fs, numel(eeg_channel_names), num_recording_samples);
if ~isempty(selected_emg_name)
    fprintf('Optional diagnostic EMG reference: %s\n', selected_emg_name);
else
    fprintf('Optional diagnostic EMG reference: none\n');
end

%% Resolve [start, end) interval
[resolved_start_sec, resolved_end_sec, start_sample, end_sample] = ...
    resolve_time_interval(start_time_sec, end_time_sec, ...
    recording_duration_sec, Fs, num_recording_samples);
selected_indices = start_sample:end_sample;
original_segment = select_structure_samples(eeg_data, selected_indices);
auxiliary_segment = select_structure_samples(auxiliary_data, selected_indices);
if isempty(selected_emg_signal)
    emg_segment = zeros(0, 1);
else
    emg_segment = double(selected_emg_signal(selected_indices));
end

fprintf('Output interval: %.6g-%.6g s (%d samples)\n', ...
    resolved_start_sec, resolved_end_sec, numel(selected_indices));
fprintf(['Windowed MTLCCA: %.6g s windows, %.0f%% overlap, ' ...
    '%d lagged EEG datasets spanning 0-%.1f ms\n'], window_duration_sec, ...
    100 * overlap_fraction, num_lag_datasets, 1000 * max_lag_duration_sec);
fprintf(['Automatic soft attenuation: at most %.0f%% candidates/window; ' ...
    'minimum automatic gain %.2f\n'], ...
    100 * max_candidate_fraction, minimum_component_gain);

%% Perform windowed multi-time-lag CCA
cca_results = get_eeg_windowed_mtlcca_artifact_removal( ...
    original_segment, Fs, ...
    WindowDurationSec=window_duration_sec, ...
    OverlapFraction=overlap_fraction, ...
    LagSamples=lag_samples, ...
    NumLagDatasets=num_lag_datasets, ...
    MaxLagDurationSec=max_lag_duration_sec, ...
    CovarianceRegularization=covariance_regularization, ...
    MCCAMaxSweeps=mcca_max_sweeps, ...
    MCCATolerance=mcca_tolerance, ...
    MCCAInitializations=mcca_initializations, ...
    MaxCandidateFraction=max_candidate_fraction, ...
    MinComponentGain=minimum_component_gain, ...
    EMGReference=emg_segment, ...
    ComponentGainsOverride=component_gains_override, ...
    ApplyPreservationGuard=apply_preservation_guard, ...
    MinLowBandCorrelation=minimum_low_band_correlation, ...
    MaxLowBandPowerChangeDb=maximum_low_band_power_change_db, ...
    MinSpatialCovarianceCorrelation=minimum_spatial_covariance_correlation, ...
    RobustClipMAD=robust_clip_mad);
cleaned_segment = cca_results.cleaned_eeg;
artifact_segment = cca_results.artifact_estimate;

num_attenuated = sum(cellfun(@(gain) nnz(gain < 1 - 1e-12), ...
    cca_results.component_gains));
num_skipped = sum(strlength(cca_results.skipped_reason) > 0);
fprintf('MTLCCA windows: %d; attenuated component instances: %d; skipped windows: %d\n', ...
    numel(cca_results.window_start_samples), num_attenuated, num_skipped);

%% Assemble a self-describing result structure
EEG_MTLCCA_Result = struct();
EEG_MTLCCA_Result.input_file = char(input_file);
EEG_MTLCCA_Result.Fs = Fs;
EEG_MTLCCA_Result.start_time_sec = resolved_start_sec;
EEG_MTLCCA_Result.end_time_sec = resolved_end_sec;
EEG_MTLCCA_Result.time_sec = ((start_sample:end_sample).' - 1) / Fs;
EEG_MTLCCA_Result.channel_names = eeg_channel_names;
EEG_MTLCCA_Result.eeg_original = original_segment;
EEG_MTLCCA_Result.eeg_clean = cleaned_segment;
EEG_MTLCCA_Result.artifact_estimate = artifact_segment;
EEG_MTLCCA_Result.emg_reference_name = selected_emg_name;
EEG_MTLCCA_Result.emg_reference = emg_segment;
EEG_MTLCCA_Result.auxiliary_signals = auxiliary_segment;
EEG_MTLCCA_Result.source_channel_labels = source_metadata.raw_labels;
EEG_MTLCCA_Result.source_channel_names = source_metadata.channel_names;
EEG_MTLCCA_Result.window_start_samples = cca_results.window_start_samples;
EEG_MTLCCA_Result.window_end_samples = cca_results.window_end_samples;
EEG_MTLCCA_Result.window_start_sec = ...
    resolved_start_sec + cca_results.window_start_sec;
EEG_MTLCCA_Result.window_end_sec = ...
    resolved_start_sec + cca_results.window_end_sec;
diagnostic_fields = {'lag_samples','lag_sec','between_lag_scores', ...
    'component_features','component_gains','component_unmixing', ...
    'component_mixing','lag_unmixing','effective_eeg_rank', ...
    'mcca_objective','mcca_sweeps','reconstruction_error', ...
    'guard_strength','preservation_before_guard', ...
    'preservation_after_guard','skipped_reason','metrics','settings','notes'};
for field_index = 1:numel(diagnostic_fields)
    field_name = diagnostic_fields{field_index};
    EEG_MTLCCA_Result.(field_name) = cca_results.(field_name);
end

%% Save numerical results, diagnostic figure, and optional EDF
[~, input_stem] = fileparts(input_file);
analysis_tag = sprintf('%s_%g-%gs', input_stem, ...
    resolved_start_sec, resolved_end_sec);
analysis_tag = regexprep(analysis_tag, '[^A-Za-z0-9_-]', '_');
output_folder = fullfile(script_dir, ...
    'EEG_MTLCCA_ArtifactRemoval_Results', analysis_tag);

write_edf_output = save_edf_result && strcmpi(source_metadata.format, 'EDF');
if save_edf_result && ~write_edf_output
    fprintf('EDF export skipped because the selected input is not EDF.\n');
end
if save_results || plot_results || write_edf_output
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
end

if write_edf_output
    proposed_edf_file = fullfile(output_folder, ...
        sprintf('%s_Windowed-MTLCCA-cleaned.edf', input_stem));
    output_edf_file = make_available_filename(proposed_edf_file);
    write_cleaned_edf(output_edf_file, source_metadata, all_signals, ...
        cleaned_segment, start_sample, end_sample, resolved_start_sec, ...
        resolved_end_sec, eeg_channel_names);
    EEG_MTLCCA_Result.output_edf_file = output_edf_file;
    fprintf('Cleaned EDF saved to: %s\n', output_edf_file);
else
    EEG_MTLCCA_Result.output_edf_file = '';
end

if save_results
    output_mat_file = fullfile(output_folder, ...
        'EEG_Windowed_MTLCCA_ArtifactRemoval_Result.mat');
    save(output_mat_file, 'EEG_MTLCCA_Result', '-v7.3');
    write_component_diagnostics_csv(EEG_MTLCCA_Result, fullfile(output_folder, ...
        'EEG_Windowed_MTLCCA_Component_Diagnostics.csv'));
    write_window_diagnostics_csv(EEG_MTLCCA_Result, fullfile(output_folder, ...
        'EEG_Windowed_MTLCCA_Window_Diagnostics.csv'));
    fprintf('Numerical results saved to: %s\n', output_mat_file);
end

if plot_results
    output_figure_file = fullfile(output_folder, ...
        'EEG_Windowed_MTLCCA_ArtifactRemoval_Dashboard.png');
    export_diagnostic_figure(EEG_MTLCCA_Result, ...
        output_figure_file, plot_resolution_dpi);
    fprintf('Diagnostic figure saved to: %s\n', output_figure_file);
end

fprintf(['Windowed MTLCCA artifact removal complete. ' ...
    'EEG_MTLCCA_Result remains in the MATLAB workspace.\n']);

%% Local loading, validation, and export helpers
function [signals, Fs, source_metadata] = ...
        load_recording(input_file, mat_channel_names)
    [~, ~, extension] = fileparts(input_file);
    switch lower(extension)
        case '.edf'
            [signals, Fs, source_metadata] = load_edf_recording(input_file);
        case '.mat'
            [signals, Fs] = load_mat_recording(input_file, mat_channel_names);
            source_metadata = struct();
            source_metadata.format = 'MAT';
            source_metadata.raw_labels = fieldnames(signals);
            source_metadata.channel_names = fieldnames(signals);
            source_metadata.info = [];
        otherwise
            error('EEG_BSS_CCA_ArtifactRemoval_runner:UnsupportedFileType', ...
                'Unsupported extension "%s". Use .edf or .mat.', extension);
    end
    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'Fs');
    Fs = double(Fs);
end

function [signals, Fs, source_metadata] = load_edf_recording(input_file)
    if exist('edfinfo', 'file') ~= 2 || exist('edfread', 'file') ~= 2
        error('EEG_BSS_CCA_ArtifactRemoval_runner:MissingEDFSupport', ...
            'MATLAB edfinfo and edfread functions are required for EDF input.');
    end

    info = edfinfo(input_file);
    record_duration_sec = seconds(info.DataRecordDuration);
    if ~isscalar(record_duration_sec) || ~isfinite(record_duration_sec) || ...
            record_duration_sec <= 0
        error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFRecordDuration', ...
            'EDF data-record duration must be a positive scalar.');
    end

    sample_rates = double(info.NumSamples(:)) / record_duration_sec;
    raw_labels = cellstr(string(info.SignalLabels(:)));
    if numel(raw_labels) ~= numel(sample_rates)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFMetadataMismatch', ...
            'EDF signal-label and sample-rate counts do not match.');
    end

    channel_names = make_unique_edf_channel_names(raw_labels);
    table_names = matlab.lang.makeValidName(raw_labels, ...
        'ReplacementStyle', 'underscore');
    [~, ~, table_name_groups] = unique(lower(string(table_names)));
    table_name_counts = accumarray(table_name_groups, 1);
    colliding = table_name_counts(table_name_groups) > 1;
    retained = ~cellfun(@is_annotation_name, raw_labels);
    signals = struct();
    retained_rates = [];

    % Read noncolliding signals in one call. Signals whose EDF labels map to
    % duplicate table names are read separately so edfread never constructs
    % a table with duplicate variables.
    bulk_indices = find(retained & ~colliding);
    if ~isempty(bulk_indices)
        bulk_table = edfread(input_file, ...
            'SelectedSignals', raw_labels(bulk_indices));
        bulk_variable_names = bulk_table.Properties.VariableNames;
        for bulk_index = 1:numel(bulk_indices)
            channel_index = bulk_indices(bulk_index);
            expected_name = table_names{channel_index};
            table_match = find(strcmpi(bulk_variable_names, expected_name), 1);
            if isempty(table_match)
                error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFSignalNotReturned', ...
                    'edfread did not return requested signal "%s".', ...
                    raw_labels{channel_index});
            end
            signal = flatten_edf_column( ...
                bulk_table.(bulk_variable_names{table_match}), ...
                raw_labels{channel_index});
            signals.(channel_names{channel_index}) = signal;
            retained_rates(end + 1) = sample_rates(channel_index); %#ok<AGROW>
        end
    end

    collision_indices = find(retained & colliding);
    for collision_index = 1:numel(collision_indices)
        channel_index = collision_indices(collision_index);
        single_table = edfread(input_file, ...
            'SelectedSignals', raw_labels{channel_index});
        if width(single_table) ~= 1
            error('EEG_BSS_CCA_ArtifactRemoval_runner:AmbiguousRawEDFLabel', ...
                ['EDF label "%s" selected %d signals. Duplicate raw EDF ' ...
                 'labels cannot be distinguished safely.'], ...
                raw_labels{channel_index}, width(single_table));
        end
        signal = flatten_edf_column( ...
            single_table.(single_table.Properties.VariableNames{1}), ...
            raw_labels{channel_index});
        signals.(channel_names{channel_index}) = signal;
        retained_rates(end + 1) = sample_rates(channel_index); %#ok<AGROW>
    end

    if isempty(retained_rates)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:NoSignalChannels', ...
            'No numeric signal channels were found in the EDF file.');
    end
    Fs = retained_rates(1);
    rate_tolerance = max(1, abs(Fs)) * 1e-10;
    if any(abs(retained_rates - Fs) > rate_tolerance)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:MixedSamplingRates', ...
            ['EDF EEG and auxiliary channels have different sampling rates. ' ...
             'Resample them to a common rate before adaptive cancellation.']);
    end

    source_metadata = struct();
    source_metadata.format = 'EDF';
    source_metadata.raw_labels = raw_labels;
    source_metadata.channel_names = channel_names;
    source_metadata.info = info;
end

function channel_names = make_unique_edf_channel_names(raw_labels)
    channel_names = cell(size(raw_labels));
    for channel_index = 1:numel(raw_labels)
        label = strtrim(raw_labels{channel_index});
        if strcmp(label, '+')
            label = 'Aux_plus';
        elseif strcmp(label, '-')
            label = 'Aux_minus';
        elseif ~isempty(regexp(label, '\s\+$', 'once'))
            label = [regexprep(label, '\s\+$', '') '_plus'];
        elseif ~isempty(regexp(label, '\s-$', 'once'))
            label = [regexprep(label, '\s-$', '') '_minus'];
        end
        channel_names{channel_index} = matlab.lang.makeValidName(label, ...
            'ReplacementStyle', 'underscore');
    end
    channel_names = matlab.lang.makeUniqueStrings( ...
        channel_names, {}, namelengthmax);
end

function signal = flatten_edf_column(column, channel_name)
    if iscell(column)
        if isempty(column) || any(~cellfun(@isnumeric, column))
            error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
                'EDF channel "%s" does not contain numeric records.', channel_name);
        end
        blocks = cellfun(@(value) double(value(:)), column, ...
            'UniformOutput', false);
        signal = vertcat(blocks{:});
    elseif isnumeric(column)
        signal = double(column(:));
    else
        error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
            'EDF channel "%s" is not numeric.', channel_name);
    end

    if isempty(signal) || ~isreal(signal) || any(~isfinite(signal))
        error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
            'EDF channel "%s" is empty, complex, or nonfinite.', channel_name);
    end
end

function [signals, Fs] = load_mat_recording(input_file, mat_channel_names)
    loaded = load(input_file);
    Fs = read_mat_sample_rate(loaded);
    signals = struct();

    if isfield(loaded, 'data') && isnumeric(loaded.data) && ...
            isreal(loaded.data) && ~isempty(mat_channel_names)
        labels = cellfun(@char, mat_channel_names(:), 'UniformOutput', false);
        matrix = double(loaded.data);
        if any(~isfinite(matrix(:)))
            error('EEG_BSS_CCA_ArtifactRemoval_runner:NonfiniteMATMatrix', ...
                'MAT variable data contains NaN or Inf values.');
        end

        if size(matrix, 2) == numel(labels) && ...
                size(matrix, 1) ~= numel(labels)
            % Already samples-by-channels.
        elseif size(matrix, 1) == numel(labels) && ...
                size(matrix, 2) ~= numel(labels)
            matrix = matrix.';
        elseif size(matrix, 1) == numel(labels) && ...
                size(matrix, 2) == numel(labels)
            error('EEG_BSS_CCA_ArtifactRemoval_runner:AmbiguousMATMatrix', ...
                ['Both dimensions of MAT variable data match mat_channel_names. ' ...
                 'Provide named vector fields instead.']);
        else
            matrix = [];
        end

        if ~isempty(matrix)
            valid_labels = make_unique_edf_channel_names(labels);
            for channel_index = 1:numel(valid_labels)
                signals.(valid_labels{channel_index}) = matrix(:, channel_index);
            end
            return;
        end
    end

    all_fields = fieldnames(loaded);
    excluded_fields = {'Fs', 'fs', 'data'};
    for field_index = 1:numel(all_fields)
        field_name = all_fields{field_index};
        if any(strcmpi(field_name, excluded_fields)) || ...
                contains(field_name, '_segment', 'IgnoreCase', true) || ...
                is_annotation_name(field_name)
            continue;
        end
        value = loaded.(field_name);
        if isnumeric(value) && isreal(value) && isvector(value) && ~isempty(value)
            if any(~isfinite(value(:)))
                error('EEG_BSS_CCA_ArtifactRemoval_runner:NonfiniteMATChannel', ...
                    'MAT channel "%s" contains NaN or Inf values.', field_name);
            end
            signals.(field_name) = double(value(:));
        end
    end

    if isempty(fieldnames(signals))
        error('EEG_BSS_CCA_ArtifactRemoval_runner:NoSignalChannels', ...
            ['No labeled numeric channel vectors were found. For matrix MAT ' ...
             'input, configure mat_channel_names to match variable data.']);
    end
end

function Fs = read_mat_sample_rate(loaded)
    if isfield(loaded, 'Fs')
        Fs = loaded.Fs;
        field_name = 'Fs';
    elseif isfield(loaded, 'fs')
        Fs = loaded.fs;
        field_name = 'fs';
    else
        error('EEG_BSS_CCA_ArtifactRemoval_runner:MissingSamplingRate', ...
            'MAT input must contain a scalar sampling rate named Fs or fs.');
    end
    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, field_name);
end

function [eeg_data, emg_data, auxiliary_data] = ...
        split_recording_channels(all_signals)
    eeg_data = struct();
    emg_data = struct();
    auxiliary_data = struct();
    channel_names = fieldnames(all_signals);

    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        if is_emg_name(channel_name)
            emg_data.(channel_name) = all_signals.(channel_name);
        elseif is_ekg_name(channel_name) || is_unnamed_auxiliary_name(channel_name)
            auxiliary_data.(channel_name) = all_signals.(channel_name);
        elseif ~is_annotation_name(channel_name)
            eeg_data.(channel_name) = all_signals.(channel_name);
        end
    end

    if isempty(fieldnames(eeg_data))
        error('EEG_BSS_CCA_ArtifactRemoval_runner:NoEEGChannels', ...
            'No EEG channels were found after excluding EMG, EKG/ECG, and annotations.');
    end
end

function tf = is_emg_name(channel_name)
    tf = contains(channel_name, 'emg', 'IgnoreCase', true);
end

function tf = is_ekg_name(channel_name)
    tf = contains(channel_name, 'ekg', 'IgnoreCase', true) || ...
        contains(channel_name, 'ecg', 'IgnoreCase', true);
end

function tf = is_unnamed_auxiliary_name(channel_name)
    tf = startsWith(channel_name, 'Aux_', 'IgnoreCase', true);
end

function tf = is_annotation_name(channel_name)
    tf = contains(channel_name, 'annotation', 'IgnoreCase', true);
end

function ordered = order_eeg_channels(eeg_data)
    preferred_order = {'Fp1_AVG', 'Fp2_AVG', 'F7_AVG', 'F3_AVG', ...
        'FZ_AVG', 'F4_AVG', 'F8_AVG', 'T3_AVG', 'C3_AVG', 'CZ_AVG', ...
        'C4_AVG', 'T4_AVG', 'T5_AVG', 'P3_AVG', 'PZ_AVG', 'P4_AVG', ...
        'T6_AVG', 'O1_AVG', 'O2_AVG'};
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
    ordered = orderfields(eeg_data, ordered_names);
end

function [num_samples, channel_names] = validate_signal_structure(data, label)
    channel_names = fieldnames(data);
    if isempty(channel_names)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EmptySignalStructure', ...
            '%s structure contains no channels.', label);
    end
    num_samples = numel(data.(channel_names{1}));
    if num_samples < 2
        error('EEG_BSS_CCA_ArtifactRemoval_runner:NotEnoughSamples', ...
            '%s channels must contain at least two samples.', label);
    end
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = data.(channel_name);
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidSignalChannel', ...
                '%s channel "%s" must be a finite real numeric vector.', ...
                label, channel_name);
        end
        if numel(signal) ~= num_samples
            error('EEG_BSS_CCA_ArtifactRemoval_runner:UnequalChannelLengths', ...
                '%s channel "%s" has %d samples; expected %d.', ...
                label, channel_name, numel(signal), num_samples);
        end
    end
end

function [selected_name, selected_signal] = ...
        select_emg_reference(emg_data, requested_name)
    channel_names = fieldnames(emg_data);
    if isempty(requested_name)
        if isscalar(channel_names)
            selected_name = channel_names{1};
            selected_signal = double(emg_data.(selected_name)(:));
            return;
        end

        [plus_index, minus_index, pair_base] = ...
            find_bipolar_emg_pair(channel_names);
        if ~isempty(plus_index)
            selected_name = sprintf('%s bipolar (+ minus -)', pair_base);
            selected_signal = double(emg_data.(channel_names{plus_index})(:)) - ...
                double(emg_data.(channel_names{minus_index})(:));
            return;
        end

        error('EEG_BSS_CCA_ArtifactRemoval_runner:AmbiguousEMGReference', ...
            ['Found %d EMG channels and no unique plus/minus pair. Set ' ...
             'emg_reference_name to exactly one of: %s'], ...
            numel(channel_names), strjoin(channel_names, ', '));
    end

    match = find(strcmpi(channel_names, char(requested_name)));
    if isempty(match)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EMGReferenceNotFound', ...
            'Requested EMG channel "%s" was not found. Available: %s', ...
            char(requested_name), strjoin(channel_names, ', '));
    end
    if numel(match) > 1
        error('EEG_BSS_CCA_ArtifactRemoval_runner:AmbiguousEMGReference', ...
            'Requested EMG name "%s" matched more than one channel.', ...
            char(requested_name));
    end
    selected_name = channel_names{match};
    selected_signal = double(emg_data.(selected_name)(:));
end

function [plus_index, minus_index, pair_base] = ...
        find_bipolar_emg_pair(channel_names)
    plus_index = [];
    minus_index = [];
    pair_base = '';
    plus_candidates = find(endsWith(channel_names, '_plus', ...
        'IgnoreCase', true));
    minus_candidates = find(endsWith(channel_names, '_minus', ...
        'IgnoreCase', true));
    matches = zeros(0, 2);
    bases = cell(0, 1);

    for plus_candidate = plus_candidates(:).'
        plus_base = regexprep(channel_names{plus_candidate}, ...
            '_plus$', '', 'ignorecase');
        for minus_candidate = minus_candidates(:).'
            minus_base = regexprep(channel_names{minus_candidate}, ...
                '_minus$', '', 'ignorecase');
            if strcmpi(plus_base, minus_base)
                matches(end + 1, :) = ...
                    [plus_candidate minus_candidate]; %#ok<AGROW>
                bases{end + 1, 1} = plus_base; %#ok<AGROW>
            end
        end
    end

    if size(matches, 1) == 1
        plus_index = matches(1, 1);
        minus_index = matches(1, 2);
        pair_base = bases{1};
    end
end

function [start_sec, end_sec, start_sample, end_sample] = ...
        resolve_time_interval(requested_start, requested_end, duration_sec, ...
        Fs, num_samples)
    if isempty(requested_start)
        start_sec = 0;
    else
        validateattributes(requested_start, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, ...
            'start_time_sec');
        start_sec = double(requested_start);
    end
    if isempty(requested_end)
        end_sec = duration_sec;
    else
        validateattributes(requested_end, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'positive'}, mfilename, ...
            'end_time_sec');
        end_sec = double(requested_end);
    end

    if end_sec <= start_sec
        error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidTimeInterval', ...
            'end_time_sec must be greater than start_time_sec.');
    end
    tolerance = max(1, duration_sec) * 1e-12;
    if end_sec > duration_sec + tolerance
        error('EEG_BSS_CCA_ArtifactRemoval_runner:IntervalExceedsRecording', ...
            'Requested end time %.6g s exceeds recording duration %.6g s.', ...
            end_sec, duration_sec);
    end
    end_sec = min(end_sec, duration_sec);

    start_sample = floor(start_sec * Fs) + 1;
    end_sample = min(floor(end_sec * Fs), num_samples);
    if end_sample < start_sample
        error('EEG_BSS_CCA_ArtifactRemoval_runner:IntervalTooShort', ...
            'The selected interval contains no samples.');
    end
end

function selected = select_structure_samples(data, indices)
    selected = struct();
    channel_names = fieldnames(data);
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = data.(channel_name);
        selected.(channel_name) = double(signal(indices(:)));
    end
end

function output_file = make_available_filename(proposed_file)
    output_file = proposed_file;
    if ~isfile(output_file)
        return;
    end

    [folder, stem, extension] = fileparts(proposed_file);
    suffix = 1;
    while isfile(output_file)
        output_file = fullfile(folder, sprintf('%s_%03d%s', ...
            stem, suffix, extension));
        suffix = suffix + 1;
    end
end

function write_cleaned_edf(output_file, source_metadata, all_signals, ...
        cleaned_eeg, start_sample, end_sample, start_time_sec, ...
        end_time_sec, eeg_channel_names)
    if exist('edfwrite', 'file') ~= 2 || exist('edfheader', 'file') ~= 2
        error('EEG_BSS_CCA_ArtifactRemoval_runner:MissingEDFWriteSupport', ...
            'MATLAB edfwrite and edfheader functions are required for EDF export.');
    end

    info = source_metadata.info;
    source_names = source_metadata.channel_names(:);
    raw_labels = string(source_metadata.raw_labels(:));
    num_signals = numel(source_names);
    num_output_samples = end_sample - start_sample + 1;
    signal_matrix = zeros(num_output_samples, num_signals);

    for signal_index = 1:num_signals
        channel_name = source_names{signal_index};
        if isfield(cleaned_eeg, channel_name)
            signal = cleaned_eeg.(channel_name);
        else
            source_signal = all_signals.(channel_name);
            signal = source_signal(start_sample:end_sample);
        end
        signal = double(signal(:));
        if numel(signal) ~= num_output_samples || any(~isfinite(signal))
            error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFOutputSignal', ...
                'Output signal "%s" has an invalid length or nonfinite samples.', ...
                channel_name);
        end
        signal_matrix(:, signal_index) = signal;
    end

    original_samples_per_record = double(info.NumSamples(:));
    samples_per_record = num_output_samples;
    for signal_index = 1:numel(original_samples_per_record)
        samples_per_record = gcd(samples_per_record, ...
            round(original_samples_per_record(signal_index)));
    end
    samples_per_record = max(1, samples_per_record);
    num_data_records = num_output_samples / samples_per_record;
    if num_data_records ~= round(num_data_records)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:InvalidEDFRecordPartition', ...
            'Output samples cannot be divided into complete EDF data records.');
    end
    num_data_records = round(num_data_records);
    record_duration_sec = samples_per_record / ...
        (double(info.NumSamples(1)) / seconds(info.DataRecordDuration));

    if startsWith(string(info.Reserved), 'EDF', 'IgnoreCase', true)
        file_type = 'EDF+';
    else
        file_type = 'EDF';
    end
    header = edfheader(file_type);
    header.Patient = string(info.Patient);
    header.Recording = string(info.Recording);
    [header.StartDate, header.StartTime] = shifted_edf_start_time( ...
        string(info.StartDate), string(info.StartTime), start_time_sec);
    header.NumDataRecords = num_data_records;
    header.DataRecordDuration = seconds(record_duration_sec);
    header.NumSignals = num_signals;
    header.SignalLabels = raw_labels;
    header.TransducerTypes = string(info.TransducerTypes(:));
    physical_dimensions = string(info.PhysicalDimensions(:));

    physical_min = double(info.PhysicalMin(:));
    physical_max = double(info.PhysicalMax(:));
    output_min = min(signal_matrix, [], 1).';
    output_max = max(signal_matrix, [], 1).';
    physical_min = min(physical_min, output_min);
    physical_max = max(physical_max, output_max);
    constant_range = physical_max <= physical_min;
    if any(constant_range)
        padding = max(1, abs(physical_min(constant_range))) * 1e-6;
        physical_min(constant_range) = physical_min(constant_range) - padding;
        physical_max(constant_range) = physical_max(constant_range) + padding;
    end
    digital_min = double(info.DigitalMin(:));
    digital_max = double(info.DigitalMax(:));
    % Keep each channel's original calibration whenever its cleaned samples
    % fit. Only the channel whose samples exceed that range is expanded.
    % This avoids degrading unaffected channels' EDF quantization.
    header.PhysicalDimensions = physical_dimensions;
    header.PhysicalMin = physical_min;
    header.PhysicalMax = physical_max;
    header.DigitalMin = digital_min;
    header.DigitalMax = digital_max;

    prefilter = string(info.Prefilter(:));
    for eeg_index = 1:numel(eeg_channel_names)
        source_index = find(strcmp(source_names, eeg_channel_names{eeg_index}), 1);
        if isempty(source_index)
            continue;
        end
        if strlength(prefilter(source_index)) == 0
            updated_prefilter = "EEG windowed MTLCCA";
        else
            updated_prefilter = prefilter(source_index) + "; MTLCCA";
        end
        if strlength(updated_prefilter) > 80
            updated_prefilter = extractBefore(updated_prefilter, 81);
        end
        prefilter(source_index) = strtrim(updated_prefilter);
    end
    header.Prefilter = prefilter;
    header.SignalReserved = string(info.SignalReserved(:));

    annotations = select_segment_annotations( ...
        info.Annotations, start_time_sec, end_time_sec);
    try
        if strcmp(file_type, 'EDF+') && ~isempty(annotations)
            edfwrite(output_file, header, signal_matrix, annotations, ...
                'InputSampleType', 'physical', ...
                'AnnotationsEncoding', 'UTF-8');
        else
            edfwrite(output_file, header, signal_matrix, ...
                'InputSampleType', 'physical');
        end
        validate_written_edf(output_file, raw_labels, num_output_samples, ...
            num_signals, physical_min, physical_max, digital_min, digital_max);
    catch exception
        if isfile(output_file)
            delete(output_file);
        end
        rethrow(exception);
    end
end

function [shifted_date, shifted_time] = shifted_edf_start_time( ...
        original_date, original_time, offset_sec)
    shifted_date = original_date;
    shifted_time = original_time;
    try
        start_datetime = datetime(original_date + " " + original_time, ...
            'InputFormat', 'dd.MM.yy HH.mm.ss', 'PivotYear', 1985);
        start_datetime = start_datetime + seconds(offset_sec);
        start_datetime.Format = 'dd.MM.yy';
        shifted_date = string(start_datetime);
        start_datetime.Format = 'HH.mm.ss';
        shifted_time = string(start_datetime);
    catch exception
        warning('EEG_BSS_CCA_ArtifactRemoval_runner:StartTimeNotShifted', ...
            ['Could not shift EDF start date/time by %.6g seconds (%s). ' ...
             'The original header time will be retained.'], ...
            offset_sec, exception.message);
    end
end

function selected = select_segment_annotations(annotations, start_sec, end_sec)
    selected = annotations;
    if isempty(annotations)
        return;
    end
    onsets = annotations.Properties.RowTimes;
    keep = ~ismissing(onsets) & ...
        onsets >= seconds(start_sec) & onsets < seconds(end_sec);
    selected = annotations(keep, :);
    selected.Properties.RowTimes = onsets(keep) - seconds(start_sec);

    % edfread represents point annotations without an explicit duration as
    % NaN seconds. edfwrite requires every duration to be non-missing, so
    % retain these as zero-duration point events.
    if ismember('Duration', selected.Properties.VariableNames)
        missing_duration = ismissing(selected.Duration);
        selected.Duration(missing_duration) = seconds(0);
    end
    if ismember('Annotations', selected.Properties.VariableNames)
        missing_text = ismissing(selected.Annotations);
        selected.Annotations(missing_text) = "Unlabeled annotation";
    end
end

function validate_written_edf(output_file, expected_labels, ...
        expected_samples, expected_signals, expected_physical_min, ...
        expected_physical_max, expected_digital_min, expected_digital_max)
    written_info = edfinfo(output_file);
    if written_info.NumSignals ~= expected_signals
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF contains %d signals; expected %d.', ...
            written_info.NumSignals, expected_signals);
    end
    if ~isequal(string(written_info.SignalLabels(:)), expected_labels(:))
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF signal labels do not match the source signal order.');
    end
    total_samples = double(written_info.NumSamples(:)) * ...
        double(written_info.NumDataRecords);
    if any(total_samples ~= expected_samples)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF signal lengths do not match the selected segment.');
    end

    written_resolution = (double(written_info.PhysicalMax(:)) - ...
        double(written_info.PhysicalMin(:))) ./ ...
        (double(written_info.DigitalMax(:)) - ...
        double(written_info.DigitalMin(:)));
    expected_resolution = (expected_physical_max(:) - expected_physical_min(:)) ./ ...
        (expected_digital_max(:) - expected_digital_min(:));
    if any(~isfinite(written_resolution)) || any(written_resolution <= 0)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF has an invalid per-channel physical resolution.');
    end
    % Header decimal formatting can introduce a tiny discrepancy. A >0.1%%
    % increase indicates avoidable quantization degradation and is rejected.
    if any(written_resolution > 1.001 * expected_resolution + eps)
        error('EEG_BSS_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF unnecessarily degraded at least one channel''s resolution.');
    end
end

function export_diagnostic_figure(result, output_file, resolution_dpi)
    ratios = result.metrics.ArtifactToSignalRMSRatio;
    finite_ratios = ratios;
    finite_ratios(~isfinite(finite_ratios)) = -Inf;
    [~, diagnostic_index] = max(finite_ratios);
    if isempty(diagnostic_index) || all(finite_ratios == -Inf)
        diagnostic_index = 1;
    end

    channel_name = result.channel_names{diagnostic_index};
    time = result.time_sec(:);
    original = result.eeg_original.(channel_name)(:);
    cleaned = result.eeg_clean.(channel_name)(:);
    artifact = result.artifact_estimate.(channel_name)(:);

    max_display_samples = max(2, round(20 * result.Fs));
    if numel(time) <= max_display_samples
        display_indices = 1:numel(time);
    else
        [~, peak_index] = max(abs(artifact));
        first_index = peak_index - floor(max_display_samples / 2);
        first_index = max(1, min(first_index, ...
            numel(time) - max_display_samples + 1));
        display_indices = first_index:(first_index + max_display_samples - 1);
    end

    window_strength = zeros(numel(result.component_gains), 1);
    for window_index = 1:numel(window_strength)
        gains = result.component_gains{window_index};
        if ~isempty(gains)
            window_strength(window_index) = max(1 - gains);
        end
    end
    [~, diagnostic_window] = max(window_strength);
    correlations = result.between_lag_scores{diagnostic_window};
    gains = result.component_gains{diagnostic_window};

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1550 1200]);
    close_figure = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout, 1);
    plot(ax1, time(display_indices), original(display_indices), ...
        'Color', [0.75 0.15 0.15]);
    hold(ax1, 'on');
    plot(ax1, time(display_indices), cleaned(display_indices), ...
        'Color', [0.1 0.35 0.75]);
    ylabel(ax1, 'EEG');
    title(ax1, sprintf( ...
        '%s: most affected diagnostic channel (not the best result)', ...
        channel_name), 'Interpreter', 'none');
    legend(ax1, {'Original', 'Cleaned'}, 'Location', 'best');
    grid(ax1, 'on');

    ax2 = nexttile(layout, 2);
    plot(ax2, time(display_indices), artifact(display_indices), ...
        'Color', [0.55 0.2 0.65]);
    ylabel(ax2, 'Removed');
    title(ax2, 'MTLCCA artifact estimate = original - cleaned');
    grid(ax2, 'on');

    ax3 = nexttile(layout, 3);
    if isempty(correlations)
        text(ax3, 0.5, 0.5, 'MTLCCA skipped in the diagnostic window', ...
            'HorizontalAlignment', 'center');
        axis(ax3, 'off');
    else
        component_numbers = (1:numel(correlations)).';
        colors = [0.85 - 0.65 * gains, 0.25 + 0.35 * gains, ...
            0.15 + 0.60 * gains];
        bars = bar(ax3, component_numbers, correlations, 'FaceColor', 'flat');
        bars.CData = colors;
        ylim(ax3, [0 1]);
        ylabel(ax3, 'Between-lag correlation');
        xlabel(ax3, 'MTLCCA component');
        attenuated = find(gains < 1 - 1e-12).';
        title(ax3, sprintf('Window %.3f-%.3f s; attenuated: %s; guard %.3f', ...
            result.window_start_sec(diagnostic_window), ...
            result.window_end_sec(diagnostic_window), ...
            mat2str(attenuated), result.guard_strength(diagnostic_window)), ...
            'Interpreter', 'none');
        grid(ax3, 'on');
    end

    ax4 = nexttile(layout, 4);
    bar(ax4, [result.metrics.HighFrequencyPowerChangeDb, ...
        result.metrics.LowBandPowerChangeDb]);
    ax4.XTick = 1:numel(result.channel_names);
    ax4.XTickLabel = result.channel_names;
    ax4.XTickLabelRotation = 45;
    ylabel(ax4, 'power change (dB)');
    xlabel(ax4, 'EEG channel');
    title(ax4, 'Spectral change after conservative attenuation');
    legend(ax4, {'30-100 Hz (58-62 excluded)', '0.5-25 Hz'}, ...
        'Location', 'best');
    grid(ax4, 'on');

    title(layout, sprintf([ ...
        'EEG windowed MTLCCA | %.3f-%.3f s | %.3g s windows | ' ...
        '%.0f%% overlap | %d lags (%d-%d samples) | minimum gain %.2f'], ...
        result.start_time_sec, result.end_time_sec, ...
        result.settings.window_duration_sec, ...
        100 * result.settings.overlap_fraction, ...
        numel(result.settings.lag_samples), min(result.settings.lag_samples), ...
        max(result.settings.lag_samples), result.settings.minimum_component_gain), ...
        'Interpreter', 'none');
    exportgraphics(fig, output_file, 'Resolution', resolution_dpi);
end

function write_component_diagnostics_csv(result, output_file)
    rows = table();
    for window_index = 1:numel(result.component_features)
        features = result.component_features{window_index};
        if isempty(features), continue; end
        count = height(features);
        window = repmat(window_index, count, 1);
        start_sec = repmat(result.window_start_sec(window_index), count, 1);
        end_sec = repmat(result.window_end_sec(window_index), count, 1);
        gain = result.component_gains{window_index}(:);
        guard_strength = repmat(result.guard_strength(window_index), count, 1);
        block = [table(window, start_sec, end_sec), features, ...
            table(gain, guard_strength)];
        rows = [rows; block]; %#ok<AGROW>
    end
    if isempty(rows)
        rows = table('Size', [0 11], ...
            'VariableTypes', repmat({'double'}, 1, 11), ...
            'VariableNames', {'window','start_sec','end_sec','Component', ...
            'Predictability','HighFrequencyFraction','SpectralSlope', ...
            'EMGCoherence','PeripheralEnergy','gain','guard_strength'});
    end
    writetable(rows, output_file);
end

function write_window_diagnostics_csv(result, output_file)
    count = numel(result.window_start_sec);
    window = (1:count).';
    start_sec = result.window_start_sec(:);
    end_sec = result.window_end_sec(:);
    effective_rank = result.effective_eeg_rank(:);
    attenuated_components = zeros(count, 1);
    minimum_gain = ones(count, 1);
    median_gain = ones(count, 1);
    for k = 1:count
        gains = result.component_gains{k};
        if isempty(gains), continue; end
        attenuated_components(k) = nnz(gains < 1 - 1e-12);
        minimum_gain(k) = min(gains);
        median_gain(k) = median(gains);
    end
    guard_strength = result.guard_strength(:);
    low_band_correlation_before = arrayfun(@(x) ...
        x.low_band_correlation, result.preservation_before_guard(:));
    low_band_power_db_before = arrayfun(@(x) ...
        x.low_band_power_change_db, result.preservation_before_guard(:));
    spatial_correlation_before = arrayfun(@(x) ...
        x.spatial_covariance_correlation, result.preservation_before_guard(:));
    low_band_correlation_after = arrayfun(@(x) ...
        x.low_band_correlation, result.preservation_after_guard(:));
    low_band_power_db_after = arrayfun(@(x) ...
        x.low_band_power_change_db, result.preservation_after_guard(:));
    spatial_correlation_after = arrayfun(@(x) ...
        x.spatial_covariance_correlation, result.preservation_after_guard(:));
    skipped_reason = result.skipped_reason(:);
    diagnostics = table(window, start_sec, end_sec, effective_rank, ...
        attenuated_components, minimum_gain, median_gain, guard_strength, ...
        low_band_correlation_before, low_band_power_db_before, ...
        spatial_correlation_before, low_band_correlation_after, ...
        low_band_power_db_after, spatial_correlation_after, skipped_reason);
    writetable(diagnostics, output_file);
end
