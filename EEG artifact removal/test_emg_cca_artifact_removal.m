function tests = test_emg_cca_artifact_removal
%TEST_EMG_CCA_ARTIFACT_REMOVAL Tests for direct EEG-EMG CCA cleanup.
%   Run with:
%       results = runtests('test_emg_cca_artifact_removal.m');
%       assertSuccess(results);

    tests = functiontests(localfunctions);
end

function testKnownEMGArtifactIsReduced(testCase)
    rng(74, 'twister');
    Fs = 200;
    num_samples = 8000;
    time = (0:num_samples - 1).' / Fs;
    emg = filter([1 -0.45 0.2], 1, randn(num_samples, 1));
    muscle = filter([0.8 -0.35 0.12], 1, emg);

    true_matrix = [ ...
        0.7 * sin(2*pi*7*time) + 0.08*randn(num_samples, 1), ...
        0.6 * sin(2*pi*11*time + 0.4) + 0.08*randn(num_samples, 1), ...
        0.5 * sin(2*pi*5*time - 0.8) + 0.08*randn(num_samples, 1), ...
        0.4 * sin(2*pi*14*time + 1.1) + 0.08*randn(num_samples, 1)];
    spatial_pattern = [1.2 -0.8 0.55 -0.35];
    observed_matrix = true_matrix + muscle * spatial_pattern;
    observed = make_structure(observed_matrix);

    results = get_emg_cca_artifact_removal(observed, emg, Fs, ...
        WindowDurationSec=10, OverlapFraction=0.5, ...
        EMGLagDurationSec=0.08, CanonicalCorrelationThreshold=0.30, ...
        MaxRejectedComponents=1);
    cleaned_matrix = structure_to_matrix(results.cleaned_eeg, ...
        results.channel_names);

    before_error = rms_value(observed_matrix - true_matrix);
    after_error = rms_value(cleaned_matrix - true_matrix);
    verifyLessThan(testCase, after_error, 0.45 * before_error);
    verifyEqual(testCase, cleaned_matrix + ...
        structure_to_matrix(results.artifact_estimate, results.channel_names), ...
        observed_matrix, 'AbsTol', 1e-10);
    verifyEqual(testCase, results.channel_names, fieldnames(observed));
    verifyTrue(testCase, all(isfinite(cleaned_matrix), 'all'));
end

function testNoRejectedComponentsReconstructExactly(testCase)
    rng(9, 'twister');
    Fs = 100;
    eeg = make_structure(randn(2300, 3));
    emg = randn(2300, 1);

    results = get_emg_cca_artifact_removal(eeg, emg, Fs, ...
        WindowDurationSec=8, OverlapFraction=0.5, ...
        CanonicalCorrelationThreshold=1, MaxRejectedComponents=1);

    verifyEqual(testCase, results.cleaned_eeg, columnize_structure(eeg), ...
        'AbsTol', 1e-10);
    verifyTrue(testCase, all(cellfun(@isempty, ...
        results.rejected_components)));
    verifyLessThan(testCase, max(results.reconstruction_error), 1e-10);
end

function testConstantEMGSkipsEveryWindow(testCase)
    rng(10, 'twister');
    eeg = make_structure(randn(1800, 3));
    emg = ones(1800, 1);

    results = get_emg_cca_artifact_removal(eeg, emg, 100, ...
        WindowDurationSec=6);

    verifyEqual(testCase, results.cleaned_eeg, columnize_structure(eeg), ...
        'AbsTol', 1e-12);
    verifyTrue(testCase, all(strlength(results.skipped_reason) > 0));
    verifyEqual(testCase, results.emg_reference_rms, 0);
end

function testAcceptsMixedVectorOrientations(testCase)
    rng(11, 'twister');
    eeg.RowChannel = randn(1, 1500);
    eeg.ColumnChannel = randn(1500, 1);
    emg = randn(1, 1500);

    results = get_emg_cca_artifact_removal(eeg, emg, 100, ...
        WindowDurationSec=5, CanonicalCorrelationThreshold=1);

    verifySize(testCase, results.cleaned_eeg.RowChannel, [1500 1]);
    verifySize(testCase, results.cleaned_eeg.ColumnChannel, [1500 1]);
end

function testRankDeficientEEGRemainsFinite(testCase)
    rng(12, 'twister');
    base = randn(2400, 1);
    eeg.ChannelA = base;
    eeg.ChannelB = base;
    eeg.ChannelC = filter(1, [1 -0.7], randn(2400, 1));
    emg = randn(2400, 1);

    results = get_emg_cca_artifact_removal(eeg, emg, 200, ...
        WindowDurationSec=6, CanonicalCorrelationThreshold=1);

    cleaned = structure_to_matrix(results.cleaned_eeg, results.channel_names);
    verifyTrue(testCase, all(isfinite(cleaned), 'all'));
    verifyTrue(testCase, all(results.effective_eeg_rank <= 2));
    verifyEqual(testCase, results.cleaned_eeg, columnize_structure(eeg), ...
        'AbsTol', 1e-10);
end

function testNumericOverrideReplacesAutomaticSelection(testCase)
    rng(13, 'twister');
    eeg = make_structure(randn(2000, 3));
    emg = randn(2000, 1);

    automatic = get_emg_cca_artifact_removal(eeg, emg, 100, ...
        WindowDurationSec=10, CanonicalCorrelationThreshold=1);
    overridden = get_emg_cca_artifact_removal(eeg, emg, 100, ...
        WindowDurationSec=10, CanonicalCorrelationThreshold=1, ...
        RejectedComponentsOverride=1);

    verifyTrue(testCase, all(cellfun(@isempty, ...
        automatic.rejected_components)));
    verifyTrue(testCase, all(cellfun(@(value) isequal(value, 1), ...
        overridden.rejected_components)));
    removed = structure_to_matrix(overridden.artifact_estimate, ...
        overridden.channel_names);
    verifyGreaterThan(testCase, norm(removed, 'fro'), 0);
end

function testCellOverrideMustMatchWindows(testCase)
    eeg = make_structure(randn(2000, 3));
    emg = randn(2000, 1);
    verifyError(testCase, @() get_emg_cca_artifact_removal( ...
        eeg, emg, 100, WindowDurationSec=5, ...
        RejectedComponentsOverride={1}), ...
        'get_emg_cca_artifact_removal:InvalidOverride');
end

function testAutomaticRejectionHonorsCap(testCase)
    rng(15, 'twister');
    emg = randn(3000, 1);
    eeg = make_structure([emg filter([1 0.4], 1, emg) ...
        filter([1 -0.3], 1, emg) randn(3000, 1)]);
    results = get_emg_cca_artifact_removal(eeg, emg, 100, ...
        WindowDurationSec=10, CanonicalCorrelationThreshold=0, ...
        MaxRejectedComponents=1);
    verifyTrue(testCase, all(cellfun(@numel, ...
        results.rejected_components) <= 1));
end

function testRejectsInvalidInputs(testCase)
    eeg = make_structure(randn(200, 3));
    emg = randn(200, 1);
    verifyError(testCase, @() get_emg_cca_artifact_removal( ...
        eeg, emg, 100, OverlapFraction=1), ...
        'get_emg_cca_artifact_removal:InvalidOverlap');
    verifyError(testCase, @() get_emg_cca_artifact_removal( ...
        eeg, emg, 100, CanonicalCorrelationThreshold=1.1), ...
        'get_emg_cca_artifact_removal:InvalidThreshold');
    verifyError(testCase, @() get_emg_cca_artifact_removal( ...
        eeg, emg(1:end-1), 100), ...
        'get_emg_cca_artifact_removal:UnequalChannelLengths');
end

function testRejectsTooShortSignal(testCase)
    eeg = make_structure(randn(20, 3));
    emg = randn(20, 1);
    verifyError(testCase, @() get_emg_cca_artifact_removal( ...
        eeg, emg, 100, WindowDurationSec=1, EMGLagDurationSec=0.18), ...
        'get_emg_cca_artifact_removal:SignalTooShort');
end

function testStatisticsToolboxIsAvailable(testCase)
    verifyEqual(testCase, exist('canoncorr', 'file'), 2, ...
        'Statistics and Machine Learning Toolbox must provide canoncorr.');
end

function structure = make_structure(matrix)
    structure = struct();
    for channel_index = 1:size(matrix, 2)
        structure.(sprintf('Channel%d', channel_index)) = ...
            matrix(:, channel_index);
    end
end

function matrix = structure_to_matrix(structure, channel_names)
    matrix = zeros(numel(structure.(channel_names{1})), numel(channel_names));
    for channel_index = 1:numel(channel_names)
        matrix(:, channel_index) = structure.(channel_names{channel_index})(:);
    end
end

function structure = columnize_structure(structure)
    names = fieldnames(structure);
    for channel_index = 1:numel(names)
        structure.(names{channel_index}) = structure.(names{channel_index})(:);
    end
end

function value = rms_value(matrix)
    value = sqrt(mean(double(matrix(:)).^2));
end
