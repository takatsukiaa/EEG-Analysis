function tests = test_eeg_windowed_mtlcca_artifact_removal
% Regression tests for conservative windowed SSQCOR MTLCCA.
    tests = functiontests(localfunctions);
end

function testExactReconstructionWithUnitGains(testCase)
    [eeg, Fs] = synthetic_eeg(6, 1600, 11);
    result = run_fast(eeg, Fs, ComponentGainsOverride=1);
    verifyLessThan(testCase, maximum_error(eeg, result.cleaned_eeg), 1e-10);
    verifyLessThan(testCase, max(result.reconstruction_error), 1e-10);
end

function testMultisetSourceRecovery(testCase)
    rng(21, 'twister');
    Fs = 200;
    time = (0:2399)' / Fs;
    sources = [sin(2*pi*4*time), filter(1, [1 -0.97], randn(size(time))), ...
        sin(2*pi*10*time + 0.2), filter([1 -1], 1, randn(size(time)))];
    sources = normalize_columns(sources);
    mixing = [1 .2 -.1 .3; .1 1 .3 -.2; -.2 .2 1 .1; .3 -.1 .1 1];
    eeg = matrix_to_struct(sources * mixing');
    result = run_fast(eeg, Fs, ComponentGainsOverride=1);
    scores = centered_matrix(eeg) * result.component_unmixing{1};
    recovery = abs(corr(scores, sources));
    verifyGreaterThan(testCase, min(max(recovery, [], 1)), 0.70);
end

function testInjectedBroadbandMuscleReduction(testCase)
    rng(31, 'twister');
    Fs = 256;
    samples = 4096;
    time = (0:samples-1)' / Fs;
    brain = [sin(2*pi*5*time), sin(2*pi*9*time+.4), ...
        filter(1, [1 -.98], randn(samples,1))];
    brain = normalize_columns(brain);
    muscle = filter([1 -1], 1, randn(samples,1));
    muscle = muscle / std(muscle);
    brain_mix = [1 .2 .1; .3 1 -.2; -.2 .4 1; .8 -.1 .3; .1 .7 .5; -.3 .2 .9];
    muscle_map = [1.4; .9; .3; .2; 1.1; .4];
    clean_truth = brain * brain_mix';
    observed = clean_truth + 1.4 * muscle * muscle_map';
    eeg = matrix_to_struct(observed, {'Fp1','Fp2','C3','C4','T3','T4'});
    result = run_fast(eeg, Fs, EMGReference=muscle, ApplyPreservationGuard=false);
    output = struct_to_matrix(result.cleaned_eeg);
    high_before = band_power(observed-clean_truth, Fs, [30 100]);
    high_after = band_power(output-clean_truth, Fs, [30 100]);
    verifyLessThan(testCase, high_after, high_before);
    verifyGreaterThan(testCase, median(abs(diag(corr(lowpass_matrix(output,Fs), ...
        lowpass_matrix(clean_truth,Fs))))), 0.90);
    verifyGreaterThan(testCase, sum(cellfun(@(x) nnz(x < 1), ...
        result.component_gains)), 0);
end

function testDeterministicConvergenceAndOrdering(testCase)
    [eeg, Fs] = synthetic_eeg(5, 1200, 41);
    first = run_fast(eeg, Fs, ComponentGainsOverride=1);
    second = run_fast(eeg, Fs, ComponentGainsOverride=1);
    verifyEqual(testCase, first.cleaned_eeg, second.cleaned_eeg, 'AbsTol', 1e-12);
    verifyEqual(testCase, first.between_lag_scores, second.between_lag_scores, ...
        'AbsTol', 1e-12);
    for k = 1:numel(first.between_lag_scores)
        verifyGreaterThanOrEqual(testCase, ...
            -diff(first.between_lag_scores{k}), -1e-12);
    end
end

function testRankDeficientAndConstantInputs(testCase)
    rng(51, 'twister');
    Fs = 200;
    x = randn(1000,1);
    eeg = matrix_to_struct([x, 2*x, -x, ones(size(x))]);
    result = run_fast(eeg, Fs, ComponentGainsOverride=1);
    verifyLessThan(testCase, maximum_error(eeg, result.cleaned_eeg), 1e-10);
    verifyLessThan(testCase, max(result.effective_eeg_rank), 4);

    constant = matrix_to_struct(ones(1000,3));
    skipped = run_fast(constant, Fs);
    verifyLessThan(testCase, maximum_error(constant, skipped.cleaned_eeg), 1e-12);
    verifyTrue(testCase, all(strlength(skipped.skipped_reason) > 0));
end

function testAbsentAndUnusableEMG(testCase)
    [eeg, Fs] = synthetic_eeg(4, 1000, 61);
    without = run_fast(eeg, Fs, ComponentGainsOverride=1);
    zeros_emg = run_fast(eeg, Fs, EMGReference=zeros(1000,1), ...
        ComponentGainsOverride=1);
    verifyEqual(testCase, without.cleaned_eeg, zeros_emg.cleaned_eeg, ...
        'AbsTol', 1e-12);
    verifyFalse(testCase, without.settings.has_emg_reference);
    verifyFalse(testCase, zeros_emg.settings.has_emg_reference);
end

function testInvalidLagVectors(testCase)
    [eeg, Fs] = synthetic_eeg(3, 800, 71);
    verifyError(testCase, @() run_fast(eeg, Fs, LagSamples=[0 2 2]), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples');
    verifyError(testCase, @() run_fast(eeg, Fs, LagSamples=[1 2 3]), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples');
    verifyError(testCase, @() run_fast(eeg, Fs, LagSamples=[0 -1 2]), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidLagSamples');
end

function testInvalidOverrides(testCase)
    [eeg, Fs] = synthetic_eeg(4, 1000, 81);
    verifyError(testCase, @() run_fast(eeg, Fs, ComponentGainsOverride=1.1), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride');
    verifyError(testCase, @() run_fast(eeg, Fs, ComponentGainsOverride=[1 .5]), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidGainOverride');
end

function testOverlapAddCoverageAndContinuity(testCase)
    [eeg, Fs] = synthetic_eeg(5, 1555, 91);
    result = run_fast(eeg, Fs, WindowDurationSec=2.1, ...
        OverlapFraction=.75, ComponentGainsOverride=1);
    verifyEqual(testCase, result.window_start_samples(1), 1);
    verifyEqual(testCase, result.window_end_samples(end), 1555);
    verifyLessThan(testCase, maximum_error(eeg, result.cleaned_eeg), 1e-10);
end

function testPreservationGuardRollsBackAggressiveGains(testCase)
    [eeg, Fs] = synthetic_eeg(6, 1600, 101);
    result = run_fast(eeg, Fs, ComponentGainsOverride=0, ...
        ApplyPreservationGuard=true);
    verifyTrue(testCase, all(result.guard_strength < 1));
    after = result.preservation_after_guard;
    verifyGreaterThanOrEqual(testCase, min([after.low_band_correlation]), .98-1e-5);
    verifyLessThanOrEqual(testCase, max(abs([after.low_band_power_change_db])), .5+1e-4);
    verifyGreaterThanOrEqual(testCase, ...
        min([after.spatial_covariance_correlation]), .95-1e-5);
end

function testInvalidAndShortInputs(testCase)
    [eeg, Fs] = synthetic_eeg(3, 800, 111);
    bad = eeg; bad.Ch1(10) = NaN;
    verifyError(testCase, @() run_fast(bad, Fs), ...
        'get_eeg_windowed_mtlcca_artifact_removal:InvalidEEGChannel');
    unequal = eeg; unequal.Ch2(end) = [];
    verifyError(testCase, @() run_fast(unequal, Fs), ...
        'get_eeg_windowed_mtlcca_artifact_removal:UnequalChannelLengths');
    short = matrix_to_struct(randn(10,4));
    verifyError(testCase, @() run_fast(short, Fs), ...
        'get_eeg_windowed_mtlcca_artifact_removal:SignalTooShort');
end

function result = run_fast(eeg, Fs, options)
    arguments
        eeg struct
        Fs double
        options.WindowDurationSec = 3
        options.OverlapFraction = .5
        options.LagSamples = [0 1 2 4]
        options.EMGReference = []
        options.ComponentGainsOverride = []
        options.ApplyPreservationGuard = true
    end
    result = get_eeg_windowed_mtlcca_artifact_removal(eeg, Fs, ...
        WindowDurationSec=options.WindowDurationSec, ...
        OverlapFraction=options.OverlapFraction, ...
        LagSamples=options.LagSamples, ...
        MCCAMaxSweeps=60, MCCATolerance=1e-7, MCCAInitializations=3, ...
        EMGReference=options.EMGReference, ...
        ComponentGainsOverride=options.ComponentGainsOverride, ...
        ApplyPreservationGuard=options.ApplyPreservationGuard);
end

function [eeg, Fs] = synthetic_eeg(channels, samples, seed)
    rng(seed, 'twister');
    Fs = 200;
    matrix = zeros(samples, channels);
    for c = 1:channels
        matrix(:,c) = filter(1, [1 -.94+.01*c], randn(samples,1));
    end
    matrix = matrix * (eye(channels) + .15*randn(channels));
    eeg = matrix_to_struct(matrix);
end

function output = matrix_to_struct(matrix, names)
    if nargin < 2
        names = arrayfun(@(x) sprintf('Ch%d',x), 1:size(matrix,2), ...
            'UniformOutput', false);
    end
    output = struct();
    for c = 1:size(matrix,2), output.(names{c}) = matrix(:,c); end
end

function matrix = struct_to_matrix(input)
    names = fieldnames(input);
    matrix = zeros(numel(input.(names{1})), numel(names));
    for c = 1:numel(names), matrix(:,c) = input.(names{c})(:); end
end

function matrix = centered_matrix(input)
    matrix = struct_to_matrix(input);
    matrix = matrix - mean(matrix,1);
end

function value = maximum_error(first, second)
    value = max(abs(struct_to_matrix(first)-struct_to_matrix(second)), [], 'all');
end

function normalized = normalize_columns(matrix)
    normalized = (matrix-mean(matrix,1))./std(matrix,0,1);
end

function filtered = lowpass_matrix(matrix, Fs)
    [b,a] = butter(4, 20/(Fs/2), 'low');
    filtered = filtfilt(b,a,matrix);
end

function power = band_power(matrix, Fs, band)
    [spectrum, frequencies] = pwelch(matrix, hamming(round(Fs)), ...
        round(Fs/2), 2^nextpow2(Fs), Fs);
    mask = frequencies >= band(1) & frequencies <= min(band(2),Fs/2-1);
    power = mean(sum(spectrum(mask,:),1));
end
