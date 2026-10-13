%% EEMD_CCA_ArtifactRemoval_runner

% Full single-channel EEMD-CCA workflow for EMG artifact removal.
% Edit only the configuration block for a different file or interval.
%
% Each EEG channel is decomposed independently. The synchronized EMG is
% supporting evidence for component rejection; delayed EEMD components are
% still the two inputs to CCA.
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
input_file = getenv('EEMD_CCA_INPUT_FILE');

% Half-open interval [start_time_sec, end_time_sec). Set either value to []
% to use the corresponding recording boundary; set both to [] for all data.
start_time_sec = [];
end_time_sec = [];

% Leave empty to select a unique EMG channel automatically or construct the
% unique detected +/− bipolar pair. Otherwise use an exact loaded field name.
emg_reference_name = '';
require_emg_reference = true;

window_duration_sec = 10;
overlap_fraction = 0.5;
ensemble_size = 10;
noise_std_ratio = 0.1;
use_complementary_noise = true;
max_num_imf = 10;
sift_relative_tolerance = 0.2;
sift_max_iterations = 100;
delay_samples = 1;

% The defaults implement the low-autocorrelation EMG formulation. The
% original motion-artifact variant can use "above" after dataset calibration.
autocorrelation_threshold = 0.9;
selection_direction = "below";
require_emg_confirmation = true;
emg_correlation_threshold = 0.10;
emg_coherence_threshold = 0.10;
emg_max_lag_sec = 0.050;
emg_coherence_band_hz = [30 95];
max_auto_rejected_fraction = 1;
artifact_component_gain = 0;

% [] uses automatic selection. A numeric vector applies globally. A cell
% array supplies indices for each window/channel (rows=windows, columns=EEG).
rejected_components_override = [];
random_seed = 0;
store_decomposition = false;

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
        'Select an EEG recording for EEMD-CCA', initial_folder);
    if isequal(selected_file, 0)
        fprintf('File selection cancelled. No data were processed.\n');
        return;
    end
    input_file = fullfile(selected_folder, selected_file);
elseif ~isfile(input_file)
    input_file = fullfile(script_dir, input_file);
end
if ~isfile(input_file)
    error('EEMD_CCA_ArtifactRemoval_runner:FileNotFound', ...
        'Input file not found: %s', input_file);
end

[all_signals, Fs, source_metadata] = load_recording( ...
    input_file, mat_channel_names);
[eeg_data, emg_data, auxiliary_data] = split_recording_channels(all_signals);
eeg_data = order_eeg_channels(eeg_data);
[num_recording_samples, eeg_channel_names] = validate_signal_structure( ...
    eeg_data, 'EEG');

selected_emg_name = ''; %#ok<NASGU>
selected_emg_signal = zeros(0, 1); %#ok<PREALL>
if ~isempty(fieldnames(emg_data))
    try
        [selected_emg_name, selected_emg_signal] = ...
            select_emg_reference(emg_data, emg_reference_name);
        if numel(selected_emg_signal) ~= num_recording_samples
            error('EEMD_CCA_ArtifactRemoval_runner:UnequalEEGAndEMGLengths', ...
                ['EEG channels contain %d samples while the EMG ' ...
                 'reference contains %d.'], num_recording_samples, ...
                numel(selected_emg_signal));
        end
    catch exception
        if require_emg_reference
            rethrow(exception);
        else
            warning('EEMD_CCA_ArtifactRemoval_runner:EMGReferenceUnavailable', ...
                ['No unique EMG reference was selected (%s). EEMD-CCA will ' ...
                 'continue with autocorrelation-only selection.'], ...
                exception.message);
            selected_emg_name = '';
            selected_emg_signal = zeros(0, 1);
        end
    end
elseif require_emg_reference
    error('EEMD_CCA_ArtifactRemoval_runner:MissingEMGReference', ...
        ['No EMG channel was found. Set require_emg_reference=false to use ' ...
         'autocorrelation-only EEMD-CCA.']);
end

recording_duration_sec = num_recording_samples / Fs;
fprintf('Loaded %s\n', input_file);
fprintf('Sampling rate: %.6g Hz; EEG channels: %d; samples: %d\n', ...
    Fs, numel(eeg_channel_names), num_recording_samples);
if ~isempty(selected_emg_name)
    fprintf('EMG selection reference: %s\n', selected_emg_name);
else
    fprintf('EMG selection reference: none (autocorrelation only)\n');
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
fprintf(['EEMD-CCA: %.6g s windows, %.0f%% overlap, %d EMD realizations, ' ...
    'noise SD %.3g, at most %d IMFs\n'], window_duration_sec, ...
    100*overlap_fraction, ensemble_size, noise_std_ratio, max_num_imf);
fprintf(['Component rule: autocorrelation %s %.3f; EMG confirmation %s; ' ...
    'artifact gain %.3f\n'], char(selection_direction), ...
    autocorrelation_threshold, string(require_emg_confirmation), ...
    artifact_component_gain);

%% Perform windowed single-channel EEMD-CCA
cca_results = get_eemd_cca_artifact_removal( ...
    original_segment, Fs, ...
    EMGReference=emg_segment, ...
    WindowDurationSec=window_duration_sec, ...
    OverlapFraction=overlap_fraction, ...
    EnsembleSize=ensemble_size, ...
    NoiseStdRatio=noise_std_ratio, ...
    UseComplementaryNoise=use_complementary_noise, ...
    MaxNumIMF=max_num_imf, ...
    SiftRelativeTolerance=sift_relative_tolerance, ...
    SiftMaxIterations=sift_max_iterations, ...
    DelaySamples=delay_samples, ...
    AutocorrelationThreshold=autocorrelation_threshold, ...
    SelectionDirection=selection_direction, ...
    RequireEMGConfirmation=require_emg_confirmation, ...
    EMGCorrelationThreshold=emg_correlation_threshold, ...
    EMGCoherenceThreshold=emg_coherence_threshold, ...
    EMGMaxLagSec=emg_max_lag_sec, ...
    EMGCoherenceBandHz=emg_coherence_band_hz, ...
    MaxAutoRejectedFraction=max_auto_rejected_fraction, ...
    ArtifactComponentGain=artifact_component_gain, ...
    RejectedComponentsOverride=rejected_components_override, ...
    RandomSeed=random_seed, ...
    StoreDecomposition=store_decomposition);
cleaned_segment = cca_results.cleaned_eeg;
artifact_segment = cca_results.artifact_estimate;

num_attenuated = sum(cellfun(@(gain) nnz(gain < 1-1e-12), ...
    cca_results.component_gains),'all');
num_skipped = sum(strlength(cca_results.skipped_reason) > 0,'all');
fprintf(['EEMD-CCA windows: %d; EEG channels: %d; attenuated component ' ...
    'instances: %d; skipped channel-windows: %d\n'], ...
    numel(cca_results.window_start_samples), numel(eeg_channel_names), ...
    num_attenuated, num_skipped);

%% Assemble a self-describing result structure
EEMD_CCA_Result = struct();
EEMD_CCA_Result.input_file = char(input_file);
EEMD_CCA_Result.Fs = Fs;
EEMD_CCA_Result.start_time_sec = resolved_start_sec;
EEMD_CCA_Result.end_time_sec = resolved_end_sec;
EEMD_CCA_Result.time_sec = ((start_sample:end_sample).' - 1)/Fs;
EEMD_CCA_Result.channel_names = eeg_channel_names;
EEMD_CCA_Result.eeg_original = original_segment;
EEMD_CCA_Result.eeg_clean = cleaned_segment;
EEMD_CCA_Result.artifact_estimate = artifact_segment;
EEMD_CCA_Result.emg_reference_name = selected_emg_name;
EEMD_CCA_Result.emg_reference = emg_segment;
EEMD_CCA_Result.auxiliary_signals = auxiliary_segment;
EEMD_CCA_Result.source_channel_labels = source_metadata.raw_labels;
EEMD_CCA_Result.source_channel_names = source_metadata.channel_names;
EEMD_CCA_Result.window_start_samples = cca_results.window_start_samples;
EEMD_CCA_Result.window_end_samples = cca_results.window_end_samples;
EEMD_CCA_Result.window_start_sec = ...
    resolved_start_sec + cca_results.window_start_sec;
EEMD_CCA_Result.window_end_sec = ...
    resolved_start_sec + cca_results.window_end_sec;
diagnostic_fields = {'component_autocorrelations','component_features', ...
    'component_gains','rejected_components','component_unmixing', ...
    'component_mixing','eemd_trial_imf_counts','eemd_component_count', ...
    'effective_component_rank','eemd_reconstruction_error', ...
    'cca_reconstruction_error','skipped_reason','metrics','settings','notes'};
if store_decomposition
    diagnostic_fields{end+1} = 'eemd_components';
end
for field_index = 1:numel(diagnostic_fields)
    field_name = diagnostic_fields{field_index};
    EEMD_CCA_Result.(field_name) = cca_results.(field_name);
end

%% Save numerical results, diagnostic figure, and optional EDF
[~, input_stem] = fileparts(input_file);
analysis_tag = sprintf('%s_%g-%gs', input_stem, ...
    resolved_start_sec, resolved_end_sec);
analysis_tag = regexprep(analysis_tag, '[^A-Za-z0-9_-]', '_');
output_folder = fullfile(script_dir, ...
    'EEMD_CCA_ArtifactRemoval_Results', analysis_tag);

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
        sprintf('%s_EEMD-CCA-cleaned.edf', input_stem));
    output_edf_file = make_available_filename(proposed_edf_file);
    write_cleaned_edf(output_edf_file, source_metadata, all_signals, ...
        cleaned_segment, start_sample, end_sample, resolved_start_sec, ...
        resolved_end_sec, eeg_channel_names);
    EEMD_CCA_Result.output_edf_file = output_edf_file;
    fprintf('Cleaned EDF saved to: %s\n', output_edf_file);
else
    EEMD_CCA_Result.output_edf_file = '';
end

if save_results
    output_mat_file = fullfile(output_folder, ...
        'EEMD_CCA_ArtifactRemoval_Result.mat');
    save(output_mat_file, 'EEMD_CCA_Result', '-v7.3');
    write_component_diagnostics_csv(EEMD_CCA_Result, fullfile(output_folder, ...
        'EEMD_CCA_Component_Diagnostics.csv'));
    write_window_diagnostics_csv(EEMD_CCA_Result, fullfile(output_folder, ...
        'EEMD_CCA_Window_Diagnostics.csv'));
    fprintf('Numerical results saved to: %s\n', output_mat_file);
end

if plot_results
    output_figure_file = fullfile(output_folder, ...
        'EEMD_CCA_ArtifactRemoval_Dashboard.png');
    export_diagnostic_figure(EEMD_CCA_Result, ...
        output_figure_file, plot_resolution_dpi);
    fprintf('Diagnostic figure saved to: %s\n', output_figure_file);
end

fprintf(['EEMD-CCA artifact removal complete. ' ...
    'EEMD_CCA_Result remains in the MATLAB workspace.\n']);

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
            error('EEMD_CCA_ArtifactRemoval_runner:UnsupportedFileType', ...
                'Unsupported extension "%s". Use .edf or .mat.', extension);
    end
    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'Fs');
    Fs = double(Fs);
end

function [signals, Fs, source_metadata] = load_edf_recording(input_file)
    if exist('edfinfo', 'file') ~= 2 || exist('edfread', 'file') ~= 2
        error('EEMD_CCA_ArtifactRemoval_runner:MissingEDFSupport', ...
            'MATLAB edfinfo and edfread functions are required for EDF input.');
    end

    info = edfinfo(input_file);
    record_duration_sec = seconds(info.DataRecordDuration);
    if ~isscalar(record_duration_sec) || ~isfinite(record_duration_sec) || ...
            record_duration_sec <= 0
        error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFRecordDuration', ...
            'EDF data-record duration must be a positive scalar.');
    end

    sample_rates = double(info.NumSamples(:)) / record_duration_sec;
    raw_labels = cellstr(string(info.SignalLabels(:)));
    if numel(raw_labels) ~= numel(sample_rates)
        error('EEMD_CCA_ArtifactRemoval_runner:EDFMetadataMismatch', ...
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
                error('EEMD_CCA_ArtifactRemoval_runner:EDFSignalNotReturned', ...
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
            error('EEMD_CCA_ArtifactRemoval_runner:AmbiguousRawEDFLabel', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:NoSignalChannels', ...
            'No numeric signal channels were found in the EDF file.');
    end
    Fs = retained_rates(1);
    rate_tolerance = max(1, abs(Fs)) * 1e-10;
    if any(abs(retained_rates - Fs) > rate_tolerance)
        error('EEMD_CCA_ArtifactRemoval_runner:MixedSamplingRates', ...
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
            error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
                'EDF channel "%s" does not contain numeric records.', channel_name);
        end
        blocks = cellfun(@(value) double(value(:)), column, ...
            'UniformOutput', false);
        signal = vertcat(blocks{:});
    elseif isnumeric(column)
        signal = double(column(:));
    else
        error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
            'EDF channel "%s" is not numeric.', channel_name);
    end

    if isempty(signal) || ~isreal(signal) || any(~isfinite(signal))
        error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFChannel', ...
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
            error('EEMD_CCA_ArtifactRemoval_runner:NonfiniteMATMatrix', ...
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
            error('EEMD_CCA_ArtifactRemoval_runner:AmbiguousMATMatrix', ...
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
                error('EEMD_CCA_ArtifactRemoval_runner:NonfiniteMATChannel', ...
                    'MAT channel "%s" contains NaN or Inf values.', field_name);
            end
            signals.(field_name) = double(value(:));
        end
    end

    if isempty(fieldnames(signals))
        error('EEMD_CCA_ArtifactRemoval_runner:NoSignalChannels', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:MissingSamplingRate', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:NoEEGChannels', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:EmptySignalStructure', ...
            '%s structure contains no channels.', label);
    end
    num_samples = numel(data.(channel_names{1}));
    if num_samples < 2
        error('EEMD_CCA_ArtifactRemoval_runner:NotEnoughSamples', ...
            '%s channels must contain at least two samples.', label);
    end
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = data.(channel_name);
        if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || ...
                isempty(signal) || any(~isfinite(signal(:)))
            error('EEMD_CCA_ArtifactRemoval_runner:InvalidSignalChannel', ...
                '%s channel "%s" must be a finite real numeric vector.', ...
                label, channel_name);
        end
        if numel(signal) ~= num_samples
            error('EEMD_CCA_ArtifactRemoval_runner:UnequalChannelLengths', ...
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

        error('EEMD_CCA_ArtifactRemoval_runner:AmbiguousEMGReference', ...
            ['Found %d EMG channels and no unique plus/minus pair. Set ' ...
             'emg_reference_name to exactly one of: %s'], ...
            numel(channel_names), strjoin(channel_names, ', '));
    end

    match = find(strcmpi(channel_names, char(requested_name)));
    if isempty(match)
        error('EEMD_CCA_ArtifactRemoval_runner:EMGReferenceNotFound', ...
            'Requested EMG channel "%s" was not found. Available: %s', ...
            char(requested_name), strjoin(channel_names, ', '));
    end
    if numel(match) > 1
        error('EEMD_CCA_ArtifactRemoval_runner:AmbiguousEMGReference', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:InvalidTimeInterval', ...
            'end_time_sec must be greater than start_time_sec.');
    end
    tolerance = max(1, duration_sec) * 1e-12;
    if end_sec > duration_sec + tolerance
        error('EEMD_CCA_ArtifactRemoval_runner:IntervalExceedsRecording', ...
            'Requested end time %.6g s exceeds recording duration %.6g s.', ...
            end_sec, duration_sec);
    end
    end_sec = min(end_sec, duration_sec);

    start_sample = floor(start_sec * Fs) + 1;
    end_sample = min(floor(end_sec * Fs), num_samples);
    if end_sample < start_sample
        error('EEMD_CCA_ArtifactRemoval_runner:IntervalTooShort', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:MissingEDFWriteSupport', ...
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
            error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFOutputSignal', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:InvalidEDFRecordPartition', ...
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
            updated_prefilter = "EEG EEMD-CCA";
        else
            updated_prefilter = prefilter(source_index) + "; EEMD-CCA";
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
        warning('EEMD_CCA_ArtifactRemoval_runner:StartTimeNotShifted', ...
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
        error('EEMD_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF contains %d signals; expected %d.', ...
            written_info.NumSignals, expected_signals);
    end
    if ~isequal(string(written_info.SignalLabels(:)), expected_labels(:))
        error('EEMD_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF signal labels do not match the source signal order.');
    end
    total_samples = double(written_info.NumSamples(:)) * ...
        double(written_info.NumDataRecords);
    if any(total_samples ~= expected_samples)
        error('EEMD_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF signal lengths do not match the selected segment.');
    end

    written_resolution = (double(written_info.PhysicalMax(:)) - ...
        double(written_info.PhysicalMin(:))) ./ ...
        (double(written_info.DigitalMax(:)) - ...
        double(written_info.DigitalMin(:)));
    expected_resolution = (expected_physical_max(:) - expected_physical_min(:)) ./ ...
        (expected_digital_max(:) - expected_digital_min(:));
    if any(~isfinite(written_resolution)) || any(written_resolution <= 0)
        error('EEMD_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
            'Written EDF has an invalid per-channel physical resolution.');
    end
    % Header decimal formatting can introduce a tiny discrepancy. A >0.1%%
    % increase indicates avoidable quantization degradation and is rejected.
    if any(written_resolution > 1.001 * expected_resolution + eps)
        error('EEMD_CCA_ArtifactRemoval_runner:EDFVerificationFailed', ...
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

    num_windows = size(result.component_gains,1);
    window_strength = zeros(num_windows,1);
    for window_index = 1:num_windows
        gains = result.component_gains{window_index,diagnostic_index};
        if ~isempty(gains)
            window_strength(window_index) = max(1-gains);
        end
    end
    [~, diagnostic_window] = max(window_strength);
    correlations = result.component_autocorrelations{ ...
        diagnostic_window,diagnostic_index};
    gains = result.component_gains{diagnostic_window,diagnostic_index};

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
    hold(ax2,'on');
    legend_entries = {'Removed EEG'};
    if ~isempty(result.emg_reference)
        emg = result.emg_reference(:);
        emg = emg-mean(emg);
        artifact_scale = max(abs(artifact(display_indices)));
        emg_scale = max(abs(emg(display_indices)));
        if artifact_scale > 0 && emg_scale > 0
            plot(ax2,time(display_indices), ...
                emg(display_indices)*(artifact_scale/emg_scale), ...
                'Color',[0.25 0.65 0.25]);
            legend_entries{end+1} = 'EMG (scaled)';
        end
    end
    legend(ax2,legend_entries,'Location','best');
    title(ax2, 'EEMD-CCA artifact estimate and synchronized EMG');
    grid(ax2, 'on');

    ax3 = nexttile(layout, 3);
    if isempty(correlations)
        text(ax3, 0.5, 0.5, 'EEMD-CCA skipped in the diagnostic channel-window', ...
            'HorizontalAlignment', 'center');
        axis(ax3, 'off');
    else
        component_numbers = (1:numel(correlations)).';
        colors = [0.85 - 0.65 * gains, 0.25 + 0.35 * gains, ...
            0.15 + 0.60 * gains];
        bars = bar(ax3, component_numbers, correlations, 'FaceColor', 'flat');
        bars.CData = colors;
        ylim(ax3, [0 1]);
        ylabel(ax3, 'Canonical autocorrelation');
        xlabel(ax3, 'EEMD-CCA component');
        attenuated = find(gains < 1 - 1e-12).';
        hold(ax3,'on');
        yline(ax3,result.settings.autocorrelation_threshold,'k--', ...
            'selection threshold');
        title(ax3, sprintf('Window %.3f-%.3f s; attenuated: %s', ...
            result.window_start_sec(diagnostic_window), ...
            result.window_end_sec(diagnostic_window), ...
            mat2str(attenuated)), ...
            'Interpreter', 'none');
        grid(ax3, 'on');
    end

    ax4 = nexttile(layout, 4);
    bar(ax4,[result.metrics.ArtifactToSignalRMSRatio, ...
        result.metrics.OriginalToCleanedCorrelation, ...
        result.metrics.OriginalMaxAbsEMGCorrelation, ...
        result.metrics.ArtifactMaxAbsEMGCorrelation]);
    ax4.XTick = 1:numel(result.channel_names);
    ax4.XTickLabel = result.channel_names;
    ax4.XTickLabelRotation = 45;
    ylabel(ax4, 'ratio / absolute correlation');
    xlabel(ax4, 'EEG channel');
    title(ax4, 'Removal magnitude, preservation, and EMG association');
    legend(ax4, {'Removed/original RMS','Original-cleaned correlation', ...
        'Original-EMG correlation','Removed-EMG correlation'}, ...
        'Location', 'best');
    grid(ax4, 'on');

    title(layout, sprintf([ ...
        'EEG EEMD-CCA | %.3f-%.3f s | %.3g s windows | ' ...
        '%.0f%% overlap | %d ensembles | noise %.3g | %s %.3f | gain %.2f'], ...
        result.start_time_sec, result.end_time_sec, ...
        result.settings.window_duration_sec, ...
        100 * result.settings.overlap_fraction, ...
        result.settings.ensemble_size, ...
        result.settings.noise_standard_deviation_ratio, ...
        result.settings.selection_direction, ...
        result.settings.autocorrelation_threshold, ...
        result.settings.artifact_component_gain), ...
        'Interpreter', 'none');
    exportgraphics(fig, output_file, 'Resolution', resolution_dpi);
end

function write_component_diagnostics_csv(result, output_file)
    rows = table();
    for window_index = 1:size(result.component_features,1)
        for channel_index = 1:size(result.component_features,2)
            features = result.component_features{window_index,channel_index};
            if isempty(features), continue; end
            count = height(features);
            Window = repmat(window_index,count,1);
            Channel = repmat(string(result.channel_names{channel_index}),count,1);
            StartSec = repmat(result.window_start_sec(window_index),count,1);
            EndSec = repmat(result.window_end_sec(window_index),count,1);
            block = [table(Window,Channel,StartSec,EndSec),features];
            rows = [rows; block]; %#ok<AGROW>
        end
    end
    if isempty(rows)
        rows = table('Size',[0 11], ...
            'VariableTypes',{'double','string','double','double','double', ...
            'double','double','double','logical','logical','double'}, ...
            'VariableNames',{'Window','Channel','StartSec','EndSec', ...
            'Component','Autocorrelation','MaxAbsEMGCorrelation', ...
            'MuscleBandEMGCoherence','EMGSupported','Selected','Gain'});
    end
    writetable(rows, output_file);
end

function write_window_diagnostics_csv(result, output_file)
    num_windows = numel(result.window_start_sec);
    num_channels = numel(result.channel_names);
    count = num_windows*num_channels;
    Window = zeros(count,1);
    Channel = strings(count,1);
    StartSec = zeros(count,1);
    EndSec = zeros(count,1);
    EEMDComponentCount = zeros(count,1);
    EffectiveRank = zeros(count,1);
    TrialIMFCountMinimum = zeros(count,1);
    TrialIMFCountMedian = zeros(count,1);
    TrialIMFCountMaximum = zeros(count,1);
    AttenuatedComponentCount = zeros(count,1);
    MinimumGain = ones(count,1);
    RejectedComponents = strings(count,1);
    EEMDReconstructionError = zeros(count,1);
    CCAReconstructionError = zeros(count,1);
    SkippedReason = strings(count,1);
    row = 0;
    for window_index = 1:num_windows
        for channel_index = 1:num_channels
            row = row+1;
            Window(row) = window_index;
            Channel(row) = string(result.channel_names{channel_index});
            StartSec(row) = result.window_start_sec(window_index);
            EndSec(row) = result.window_end_sec(window_index);
            EEMDComponentCount(row) = ...
                result.eemd_component_count(window_index,channel_index);
            EffectiveRank(row) = ...
                result.effective_component_rank(window_index,channel_index);
            trial_counts = result.eemd_trial_imf_counts{window_index,channel_index};
            if ~isempty(trial_counts)
                TrialIMFCountMinimum(row) = min(trial_counts);
                TrialIMFCountMedian(row) = median(trial_counts);
                TrialIMFCountMaximum(row) = max(trial_counts);
            end
            gains = result.component_gains{window_index,channel_index};
            if ~isempty(gains)
                AttenuatedComponentCount(row) = nnz(gains < 1-1e-12);
                MinimumGain(row) = min(gains);
            end
            rejected = result.rejected_components{window_index,channel_index};
            RejectedComponents(row) = string(mat2str(rejected));
            EEMDReconstructionError(row) = ...
                result.eemd_reconstruction_error(window_index,channel_index);
            CCAReconstructionError(row) = ...
                result.cca_reconstruction_error(window_index,channel_index);
            SkippedReason(row) = result.skipped_reason(window_index,channel_index);
        end
    end
    diagnostics = table(Window,Channel,StartSec,EndSec,EEMDComponentCount, ...
        EffectiveRank,TrialIMFCountMinimum,TrialIMFCountMedian, ...
        TrialIMFCountMaximum,AttenuatedComponentCount,MinimumGain, ...
        RejectedComponents,EEMDReconstructionError,CCAReconstructionError, ...
        SkippedReason);
    writetable(diagnostics, output_file);
end
