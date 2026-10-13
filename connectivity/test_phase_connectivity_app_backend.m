function tests = test_phase_connectivity_app_backend
%TEST_PHASE_CONNECTIVITY_APP_BACKEND Tests for the GUI analysis backend.
    tests = functiontests(localfunctions);
end

function testLoadsNamedMatChannelsAndMetadata(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    Fs = 200; %#ok<NASGU>
    Fp1_AVG = (1:1000).'; %#ok<NASGU>
    C3_AVG = (1001:2000).'; %#ok<NASGU>
    EMG1 = zeros(1000, 1); %#ok<NASGU>
    C3_segment1 = zeros(50, 1); %#ok<NASGU>
    file = fullfile(folder, 'recording.mat');
    save(file, 'Fs', 'Fp1_AVG', 'C3_AVG', 'EMG1', 'C3_segment1');
    recording = load_phase_connectivity_recording(file);
    verifyEqual(testCase, recording.Fs, 200);
    verifyEqual(testCase, recording.channel_names, {'Fp1_AVG';'C3_AVG'});
    verifyEqual(testCase, recording.duration_sec, 5);
end

function testLoadsMatrixAndChanlocLabels(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    Fs = 100; %#ok<NASGU>
    data = randn(2, 500); %#ok<NASGU>
    chanlocs = struct('labels', {'F3','F4'}, 'theta', {-39,39}, 'radius', {.5,.5}); %#ok<NASGU>
    file = fullfile(folder, 'matrix.mat'); save(file, 'Fs', 'data', 'chanlocs');
    recording = load_phase_connectivity_recording(file);
    verifyEqual(testCase, recording.channel_names, {'F3';'F4'});
    verifyEqual(testCase, recording.num_samples, 500);
    verifyEqual(testCase, numel(recording.chanlocs), 2);
end

function testBandpassAndNotchAttenuation(testCase)
    Fs = 500; time = (0:(10*Fs-1)).'/Fs;
    signal = sin(2*pi*10*time) + sin(2*pi*50*time);
    filtered = apply_eeg_filter_chain(signal, Fs, struct( ...
        'bandpass_enabled', true, 'band_hz', [5 80], ...
        'bandpass_order', 100, 'notch_hz', 50));
    [rawPsd, frequency] = periodogram(signal, [], [], Fs);
    filteredPsd = periodogram(filtered, [], [], Fs);
    [~, at10] = min(abs(frequency - 10));
    [~, at50] = min(abs(frequency - 50));
    verifyGreaterThan(testCase, filteredPsd(at10), rawPsd(at10) * .5);
    verifyLessThan(testCase, filteredPsd(at50), rawPsd(at50) * .02);
    verifySize(testCase, filtered, size(signal));
end

function testConnectivityUsesHalfOpenIntervalAndSelectedChannels(testCase)
    recording = synthetic_recording();
    options = struct('start_time_sec', 1, 'end_time_sec', 3, ...
        'channel_names', {{'Fp1_AVG','Fp2_AVG'}}, ...
        'window_length_sec', .5, 'window_step_sec', .25, ...
        'filter_settings', struct('bandpass_enabled', true, ...
        'band_hz', [5 30], 'bandpass_order', 20, 'notch_hz', 0));
    result = run_phase_connectivity_analysis(recording, options);
    verifyEqual(testCase, result.start_sample, 201);
    verifyEqual(testCase, result.end_sample, 600);
    verifyEqual(testCase, result.channel_names, {'Fp1_AVG';'Fp2_AVG'});
    verifySize(testCase, result.plv_matrix, [2 2]);
    verifyGreaterThan(testCase, result.plv_matrix(1,2), .98);
end

function testSingleChannelPsdFindsKnownPeak(testCase)
    recording = synthetic_recording();
    result = compute_single_channel_psd(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 5, 'channel_name', 'Fp1_AVG', ...
        'max_frequency_hz', 70, 'filter_settings', struct( ...
        'bandpass_enabled', false, 'band_hz', [1 70], ...
        'bandpass_order', 20, 'notch_hz', 0)));
    [~, peak] = max(result.raw_psd);
    verifyEqual(testCase, result.frequency(peak), 10, 'AbsTol', 1);
    verifyEqual(testCase, result.frequency, result.frequency(:));
end

function testPsdTopographyBandPowerAndEpochTiming(testCase)
    recording = synthetic_recording();
    result = compute_psd_topography(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 4.9, ...
        'band_hz', [9 11], 'epoch_length_sec', 0.6, ...
        'filter_settings', struct('bandpass_enabled', false, ...
        'band_hz', [1 70], 'bandpass_order', 20, 'notch_hz', 0)));
    verifyEqual(testCase, result.num_epochs, 8);
    verifyEqual(testCase, result.unused_tail_sec, 0.1, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.epoch_start_time_sec(1), 0);
    verifyEqual(testCase, result.epoch_end_time_sec(end), 4.8, 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, result.raw_whole_band_power(1), ...
        result.raw_whole_band_power(3) * 5);
    verifyEqual(testCase, result.raw_whole_band_power, ...
        mean(result.raw_band_power, 1), 'AbsTol', 1e-12);
    verifyFalse(testCase, isfield(result, 'raw_signal'));
    verifyFalse(testCase, isfield(result, 'filtered_signal'));
end

function testPsdTopographyFilteredPowerAndValidation(testCase)
    recording = synthetic_recording();
    result = compute_psd_topography(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 5, ...
        'band_hz', [9 11], 'epoch_length_sec', 1, ...
        'filter_settings', struct('bandpass_enabled', true, ...
        'band_hz', [20 30], 'bandpass_order', 20, 'notch_hz', 0)));
    verifyLessThan(testCase, result.filtered_whole_band_power(1), ...
        result.raw_whole_band_power(1) * 0.05);
    verifyError(testCase, @() compute_psd_topography(recording, struct( ...
        'band_hz', [20 101])), 'compute_psd_topography:InvalidBand');
    verifyError(testCase, @() compute_psd_topography(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 0.2, 'epoch_length_sec', 0.5)), ...
        'compute_psd_topography:IntervalTooShort');
end

function testBuiltInLocationMappingNormalizesAvg(testCase)
    mapping = resolve_electrode_locations({'Fp1_AVG','FZ-AVG','unknown'}, []);
    verifyEqual(testCase, mapping.num_mapped, 2);
    verifyEqual(testCase, mapping.analysis_indices, [1;2]);
    verifyEqual(testCase, mapping.unmatched_channels, {'unknown'});
end

function testImportsChanlocsFromSetCompatibleMat(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    chanlocs = struct('labels', {'X1','X2'}, 'theta', {-45,45}, 'radius', {.5,.5}); %#ok<NASGU>
    file = fullfile(folder, 'montage.set'); save(file, 'chanlocs');
    mapping = resolve_electrode_locations({'X2','X1'}, file);
    verifyEqual(testCase, mapping.num_mapped, 2);
    verifyEqual(testCase, mapping.analysis_indices, [1;2]);
    verifyGreaterThan(testCase, mapping.x(1), 0);
end

function testTopoplotCountsConnections(testCase)
    mapping = resolve_electrode_locations({'Fp1','Fp2','Cz'}, []);
    fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig);
    stats = topoplot_connect(struct('chanPairs', [1 2;2 3], ...
        'connectStrength', [.6;.8], 'connectStrengthLimits', [0 1]), mapping, ax);
    verifyEqual(testCase, stats.numConnections, 2);
    verifyEqual(testCase, stats.minimumStrength, .6, 'AbsTol', eps);
end

function testPsdTopographyRendersInUIAxes(testCase)
    mapping = resolve_electrode_locations({'Fp1','Fp2','C3'}, []);
    fig = uifigure('Visible','off'); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = uiaxes(fig);
    stats = plot_psd_topography([1 2 3], mapping, ax, struct( ...
        'map_limits', [0 4], 'colorbar_label', 'Test power'));
    verifyEqual(testCase, stats.num_mapped, 3);
    verifyEqual(testCase, stats.map_limits, [0 4]);
    verifyNotEmpty(testCase, findobj(ax, 'Type', 'Surface'));
end

function testExportProducesNumericsAndFigures(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    recording = synthetic_recording();
    connectivity = run_phase_connectivity_analysis(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 2, ...
        'window_length_sec', .5, 'window_step_sec', .25, ...
        'filter_settings', struct('bandpass_enabled', true, ...
        'band_hz', [5 30], 'bandpass_order', 20, 'notch_hz', 0)));
    psd = compute_single_channel_psd(recording, struct('start_time_sec',0, ...
        'end_time_sec',2,'channel_name','Fp1_AVG'));
    mapping = resolve_electrode_locations(connectivity.channel_names, []);
    files = export_phase_connectivity_results(struct('connectivity', connectivity, ...
        'psd', psd, 'electrode_mapping', mapping), folder, struct());
    verifyGreaterThanOrEqual(testCase, numel(files), 4);
    verifyTrue(testCase, all(cellfun(@isfile, files)));
    saved = load(fullfile(folder, 'PhaseConnectivityResults.mat'), 'bundle');
    verifyFalse(testCase, isfield(saved.bundle.connectivity, 'filtered_eeg'));
    verifyFalse(testCase, isfield(saved.bundle.psd, 'raw_signal'));
end

function testExportsNarrowPsdDisplayRange(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    recording = synthetic_recording();
    psd = compute_single_channel_psd(recording, struct('start_time_sec',0, ...
        'end_time_sec',5,'channel_name','Fp1_AVG','max_frequency_hz',100));
    options = struct('export_connectivity', false, 'export_topoplot', false, ...
        'export_psd', true, 'psd_source', 'Raw', ...
        'psd_frequency_range_hz', [8 20]);
    files = export_phase_connectivity_results(struct('psd', psd), folder, options);
    verifyTrue(testCase, isfile(fullfile(folder, 'Fp1_AVG_PSD_Raw.png')));
    verifyEqual(testCase, numel(files), 2);
end

function testExportsPsdTopographyCurrentAndAllEpochs(testCase)
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    recording = synthetic_recording();
    topography = compute_psd_topography(recording, struct( ...
        'start_time_sec', 0, 'end_time_sec', 2, ...
        'band_hz', [9 11], 'epoch_length_sec', 0.5));
    mapping = resolve_electrode_locations(topography.channel_names, []);
    bundle = struct('psd_topography', topography, ...
        'psd_topography_mapping', mapping);
    options = struct('export_connectivity', false, 'export_psd', false, ...
        'export_topoplot', false, 'export_psd_topography', true, ...
        'export_all_psd_topography_epochs', true, ...
        'psd_topography_source', 'Raw', 'psd_topography_in_db', true, ...
        'psd_topography_scope', 'Epoch', 'psd_topography_epoch_index', 2);
    files = export_phase_connectivity_results(bundle, folder, options);
    verifyEqual(testCase, numel(files), 3);
    verifyTrue(testCase, all(cellfun(@isfile, files)));
    verifyTrue(testCase, any(contains(files, 'Epoch_0002')));
    verifyTrue(testCase, any(contains(files, 'AllEpochs')));
    saved = load(fullfile(folder, 'PhaseConnectivityResults.mat'), 'bundle');
    verifyTrue(testCase, isfield(saved.bundle, 'psd_topography'));
    verifyFalse(testCase, isfield(saved.bundle.psd_topography, 'raw_signal'));
end

function testAppConstructs(testCase)
    app = PhaseConnectivityApp; cleanup = onCleanup(@() delete(app)); %#ok<NASGU>
    drawnow;
    verifyTrue(testCase, isvalid(app.UIFigure));
    verifyEqual(testCase, app.UIFigure.Name, 'EEG Phase Connectivity');
    tabs = findall(app.UIFigure, 'Type', 'uitab');
    verifyTrue(testCase, any(strcmp({tabs.Title}, 'PSD Topoplot')));
end

function recording = synthetic_recording()
    Fs = 200;
    time = (0:(5*Fs-1)).'/Fs;
    data.Fp1_AVG = sin(2*pi*10*time);
    data.Fp2_AVG = sin(2*pi*10*time + pi/4);
    data.C3_AVG = sin(2*pi*13*time);
    recording = struct('input_file','synthetic.mat','file_name','synthetic.mat', ...
        'extension','.mat','data',data,'Fs',Fs, ...
        'channel_names',{{'Fp1_AVG';'Fp2_AVG';'C3_AVG'}}, ...
        'num_samples',numel(time),'duration_sec',numel(time)/Fs, ...
        'chanlocs',struct([]),'source_metadata',struct('file_name','synthetic.mat'));
end
