function tests = test_emg_adaptive_cancellation
%TEST_EMG_ADAPTIVE_CANCELLATION Synthetic tests for normalized-LMS cleanup.
%   Run with:
%       results = runtests('test_emg_adaptive_cancellation.m');
%       assertSuccess(results);

    tests = functiontests(localfunctions);
end

function testKnownEMGArtifactIsReduced(testCase)
    rng(41, 'twister');
    Fs = 200;
    num_samples = 12000;
    time = (0:num_samples - 1).' / Fs;
    emg = filter([1 -0.6 0.25], 1, randn(num_samples, 1));

    true_eeg.ChannelA = 0.4 * sin(2 * pi * 9 * time) + ...
        0.08 * randn(num_samples, 1);
    true_eeg.ChannelB = 0.3 * sin(2 * pi * 13 * time + 0.7) + ...
        0.08 * randn(num_samples, 1);
    true_eeg.ChannelC = 0.2 * sin(2 * pi * 6 * time - 0.4) + ...
        0.08 * randn(num_samples, 1);

    observed.ChannelA = true_eeg.ChannelA + filter([0.9 -0.35 0.1], 1, emg);
    observed.ChannelB = true_eeg.ChannelB + filter([-0.5 0.25], 1, emg);
    observed.ChannelC = true_eeg.ChannelC + filter([0.25 0.1 -0.08], 1, emg);

    results = get_emg_adaptive_cancellation(observed, emg, Fs, ...
        FilterLengthSec=0.08, StepSize=0.05);

    burn_in = 2001:num_samples;
    channel_names = fieldnames(observed);
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        before_error = observed.(channel_name)(burn_in) - ...
            true_eeg.(channel_name)(burn_in);
        after_error = results.cleaned_eeg.(channel_name)(burn_in) - ...
            true_eeg.(channel_name)(burn_in);
        verifyLessThan(testCase, rms_value(after_error), ...
            0.25 * rms_value(before_error));
    end

    verifyEqual(testCase, results.channel_names, channel_names);
    verifySize(testCase, results.final_weights, [16 3]);
    verifyEqual(testCase, height(results.metrics), 3);
    verifyTrue(testCase, all(isfinite(results.final_weights), 'all'));
end

function testNoContaminationRemainsZero(testCase)
    rng(7, 'twister');
    Fs = 200;
    emg = randn(1000, 1);
    eeg.ZeroChannel = zeros(1, 1000);

    results = get_emg_adaptive_cancellation(eeg, emg, Fs);

    verifyEqual(testCase, results.cleaned_eeg.ZeroChannel, zeros(1000, 1), ...
        'AbsTol', 1e-14);
    verifyEqual(testCase, results.artifact_estimate.ZeroChannel, ...
        zeros(1000, 1), 'AbsTol', 1e-14);
end

function testAcceptsMixedVectorOrientations(testCase)
    rng(12, 'twister');
    Fs = 100;
    emg = randn(1, 1000);
    eeg.RowChannel = randn(1, 1000);
    eeg.ColumnChannel = randn(1000, 1);

    results = get_emg_adaptive_cancellation(eeg, emg, Fs);

    verifySize(testCase, results.cleaned_eeg.RowChannel, [1000 1]);
    verifySize(testCase, results.cleaned_eeg.ColumnChannel, [1000 1]);
    verifyTrue(testCase, all(isfinite(results.cleaned_eeg.RowChannel)));
    verifyTrue(testCase, all(isfinite(results.cleaned_eeg.ColumnChannel)));
end

function testQuietReferenceProducesNoCancellation(testCase)
    rng(18, 'twister');
    Fs = 200;
    num_samples = 2000;
    emg = [randn(300, 1); zeros(num_samples - 300, 1)];
    eeg.ChannelA = randn(num_samples, 1);

    results = get_emg_adaptive_cancellation(eeg, emg, Fs, ...
        FilterLengthSec=0.08, StepSize=0.05);

    quiet_indices = 317:num_samples;
    verifyEqual(testCase, ...
        results.artifact_estimate.ChannelA(quiet_indices), ...
        zeros(numel(quiet_indices), 1), 'AbsTol', 1e-14);
    verifyEqual(testCase, results.cleaned_eeg.ChannelA(quiet_indices), ...
        eeg.ChannelA(quiet_indices), 'AbsTol', 1e-14);
end

function testRejectsInvalidStepSize(testCase)
    eeg.ChannelA = randn(100, 1);
    emg = randn(100, 1);
    verifyError(testCase, @() get_emg_adaptive_cancellation( ...
        eeg, emg, 100, StepSize=2), ...
        'get_emg_adaptive_cancellation:InvalidStepSize');
end

function testRejectsUnequalLengths(testCase)
    eeg.ChannelA = randn(99, 1);
    emg = randn(100, 1);
    verifyError(testCase, @() get_emg_adaptive_cancellation(eeg, emg, 100), ...
        'get_emg_adaptive_cancellation:UnequalChannelLengths');
end

function testRejectsShortSignal(testCase)
    eeg.ChannelA = randn(8, 1);
    emg = randn(8, 1);
    verifyError(testCase, @() get_emg_adaptive_cancellation( ...
        eeg, emg, 100, FilterLengthSec=0.08), ...
        'get_emg_adaptive_cancellation:SignalTooShort');
end

function testRejectsConstantEMGReference(testCase)
    eeg.ChannelA = randn(100, 1);
    emg = ones(100, 1);
    verifyError(testCase, @() get_emg_adaptive_cancellation(eeg, emg, 100), ...
        'get_emg_adaptive_cancellation:ConstantEMGReference');
end

function testDSPSystemToolboxIsAvailable(testCase)
    % The implementation also raises MissingDSPToolbox when this class is
    % absent. This environment-dependent prerequisite is checked directly.
    verifyEqual(testCase, exist('dsp.LMSFilter', 'class'), 8, ...
        'DSP System Toolbox must provide dsp.LMSFilter.');
end

function value = rms_value(signal)
    value = sqrt(mean(double(signal(:)).^2));
end
