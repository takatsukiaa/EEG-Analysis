function tests = test_eeg_bss_cca_artifact_removal
%TEST_EEG_BSS_CCA_ARTIFACT_REMOVAL Tests for delayed-EEG BSS-CCA cleanup.
%   Run with:
%       results = runtests('test_eeg_bss_cca_artifact_removal.m');
%       assertSuccess(results);

    tests = functiontests(localfunctions);
end

function testLowAutocorrelationMuscleSourceIsReduced(testCase)
    rng(81, 'twister');
    Fs = 200;
    num_samples = 8000;
    time = (0:num_samples - 1).' / Fs;

    brain_sources = [ ...
        sin(2*pi*6*time) + 0.04*randn(num_samples, 1), ...
        0.8*sin(2*pi*10*time + 0.5) + 0.04*randn(num_samples, 1), ...
        filter(1, [1 -0.97], 0.03*randn(num_samples, 1))];
    muscle_source = filter([1 -0.15], 1, randn(num_samples, 1));
    mixing = orth(randn(4));
    true_matrix = [brain_sources zeros(num_samples, 1)] * mixing.';
    observed_matrix = [brain_sources 1.8*muscle_source] * mixing.';
    observed = make_structure(observed_matrix);

    results = get_eeg_bss_cca_artifact_removal(observed, Fs, ...
        WindowDurationSec=10, OverlapFraction=0.5, DelaySamples=1, ...
        AutocorrelationThreshold=0.80, MaxAutoRejectedFraction=0.5);
    cleaned_matrix = structure_to_matrix(results.cleaned_eeg, ...
        results.channel_names);

    before_error = rms_value(observed_matrix - true_matrix);
    after_error = rms_value(cleaned_matrix - true_matrix);
    verifyLessThan(testCase, after_error, 0.35 * before_error);
    verifyEqual(testCase, cleaned_matrix + ...
        structure_to_matrix(results.artifact_estimate, results.channel_names), ...
        observed_matrix, 'AbsTol', 1e-10);
    verifyTrue(testCase, all(isfinite(cleaned_matrix), 'all'));
end

function testNoRejectedComponentsReconstructExactly(testCase)
    rng(82, 'twister');
    Fs = 100;
    num_samples = 2300;
    smooth = filter(1, [1 -0.98], randn(num_samples, 3));
    eeg = make_structure(smooth);

    results = get_eeg_bss_cca_artifact_removal(eeg, Fs, ...
        WindowDurationSec=8, OverlapFraction=0.5, ...
        AutocorrelationThreshold=0, MaxAutoRejectedFraction=0.25);

    verifyEqual(testCase, results.cleaned_eeg, columnize_structure(eeg), ...
        'AbsTol', 1e-10);
    verifyTrue(testCase, all(cellfun(@isempty, ...
        results.rejected_components)));
    verifyLessThan(testCase, max(results.reconstruction_error), 1e-10);
end

function testConstantEEGSkipsEveryWindow(testCase)
    eeg.ChannelA = ones(1800, 1);
    eeg.ChannelB = 2 * ones(1800, 1);
    eeg.ChannelC = -ones(1800, 1);

    results = get_eeg_bss_cca_artifact_removal(eeg, 100, ...
        WindowDurationSec=6);

    verifyEqual(testCase, results.cleaned_eeg, eeg, 'AbsTol', 1e-12);
    verifyTrue(testCase, all(strlength(results.skipped_reason) > 0));
end

function testAcceptsMixedVectorOrientations(testCase)
    rng(83, 'twister');
    eeg.RowChannel = randn(1, 1500);
    eeg.ColumnChannel = randn(1500, 1);

    results = get_eeg_bss_cca_artifact_removal(eeg, 100, ...
        WindowDurationSec=5, AutocorrelationThreshold=0);

    verifySize(testCase, results.cleaned_eeg.RowChannel, [1500 1]);
    verifySize(testCase, results.cleaned_eeg.ColumnChannel, [1500 1]);
end

function testRankDeficientEEGRemainsFinite(testCase)
    rng(84, 'twister');
    base = filter(1, [1 -0.9], randn(2400, 1));
    eeg.ChannelA = base;
    eeg.ChannelB = base;
    eeg.ChannelC = filter(1, [1 -0.7], randn(2400, 1));

    results = get_eeg_bss_cca_artifact_removal(eeg, 200, ...
        WindowDurationSec=6, AutocorrelationThreshold=0);

    cleaned = structure_to_matrix(results.cleaned_eeg, results.channel_names);
    verifyTrue(testCase, all(isfinite(cleaned), 'all'));
    verifyTrue(testCase, all(results.effective_eeg_rank <= 2));
    verifyEqual(testCase, results.cleaned_eeg, columnize_structure(eeg), ...
        'AbsTol', 1e-10);
end

function testNumericOverrideReplacesAutomaticSelection(testCase)
    rng(85, 'twister');
    eeg = make_structure(randn(2000, 3));

    automatic = get_eeg_bss_cca_artifact_removal(eeg, 100, ...
        WindowDurationSec=10, AutocorrelationThreshold=0);
    overridden = get_eeg_bss_cca_artifact_removal(eeg, 100, ...
        WindowDurationSec=10, AutocorrelationThreshold=0, ...
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
    verifyError(testCase, @() get_eeg_bss_cca_artifact_removal( ...
        eeg, 100, WindowDurationSec=5, ...
        RejectedComponentsOverride={1}), ...
        'get_eeg_bss_cca_artifact_removal:InvalidOverride');
end

function testAutomaticRejectionHonorsFraction(testCase)
    rng(86, 'twister');
    eeg = make_structure(randn(3000, 8));
    results = get_eeg_bss_cca_artifact_removal(eeg, 100, ...
        WindowDurationSec=10, AutocorrelationThreshold=1, ...
        MaxAutoRejectedFraction=0.25);
    verifyTrue(testCase, all(cellfun(@numel, ...
        results.rejected_components) <= 2));
end

function testRejectsInvalidInputs(testCase)
    eeg = make_structure(randn(200, 3));
    verifyError(testCase, @() get_eeg_bss_cca_artifact_removal( ...
        eeg, 100, OverlapFraction=1), ...
        'get_eeg_bss_cca_artifact_removal:InvalidOverlap');
    verifyError(testCase, @() get_eeg_bss_cca_artifact_removal( ...
        eeg, 100, AutocorrelationThreshold=1.1), ...
        'get_eeg_bss_cca_artifact_removal:InvalidThreshold');
    verifyError(testCase, @() get_eeg_bss_cca_artifact_removal( ...
        eeg, 100, MaxAutoRejectedFraction=1.1), ...
        'get_eeg_bss_cca_artifact_removal:InvalidRejectedFraction');
end

function testRejectsTooShortSignal(testCase)
    eeg = make_structure(randn(10, 3));
    verifyError(testCase, @() get_eeg_bss_cca_artifact_removal( ...
        eeg, 100, WindowDurationSec=1, DelaySamples=8), ...
        'get_eeg_bss_cca_artifact_removal:SignalTooShort');
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
