function tests = test_phase_connectivity
%TEST_PHASE_CONNECTIVITY Synthetic tests for aligned PLV/PLI/wPLI analysis.
%   Run with:
%       results = runtests('test_phase_connectivity.m');
%       assertSuccess(results);

    tests = functiontests(localfunctions);
end

function testConstantNonzeroPhaseOffset(testCase)
    Fs = 500;
    time = (0:(10 * Fs - 1)).' / Fs;
    signals.ChannelA = cos(2 * pi * 12 * time);
    signals.ChannelB = cos(2 * pi * 12 * time + pi / 4);

    results = get_phase_connectivity(signals, Fs, 1, 0.5);

    verifyGreaterThan(testCase, results.plv_matrix(1, 2), 0.99);
    verifyGreaterThan(testCase, results.pli_matrix(1, 2), 0.99);
    verifyGreaterThan(testCase, results.wpli_matrix(1, 2), 0.99);
end

function testIdenticalZeroLagSignals(testCase)
    Fs = 500;
    time = (0:(4 * Fs - 1)).' / Fs;
    signal = cos(2 * pi * 10 * time) + 0.2 * cos(2 * pi * 18 * time);
    signals.ChannelA = signal;
    signals.ChannelB = signal;

    results = get_phase_connectivity(signals, Fs, 0.5, 0.25);

    verifyEqual(testCase, results.plv_matrix(1, 2), 1, 'AbsTol', 1e-14);
    verifyEqual(testCase, results.pli_matrix(1, 2), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, results.wpli_matrix(1, 2), 0, 'AbsTol', 1e-14);
end

function testIndependentSignalsAreNearZero(testCase)
    rng(29, 'twister');
    Fs = 500;
    num_samples = 50000;
    signals.ChannelA = randn(num_samples, 1);
    signals.ChannelB = randn(num_samples, 1);

    results = get_phase_connectivity(signals, Fs, 2, 1);

    verifyLessThan(testCase, results.plv_matrix(1, 2), 0.08);
    verifyLessThan(testCase, results.pli_matrix(1, 2), 0.08);
    verifyLessThan(testCase, results.wpli_matrix(1, 2), 0.08);
end

function testMatrixContractsAndAlignment(testCase)
    Fs = 100;
    num_samples = 1000;
    time = (0:(num_samples - 1)).' / Fs;
    signals.C3 = cos(2 * pi * 7 * time);
    signals.C4 = cos(2 * pi * 7 * time + pi / 3);
    signals.Pz = cos(2 * pi * 9 * time);

    results = get_phase_connectivity(signals, Fs, 0.2, 0.1);

    verifyEqual(testCase, results.channel_names, {'C3'; 'C4'; 'Pz'});
    verifySize(testCase, results.plv_matrix, [3 3]);
    verifySize(testCase, results.pli_matrix, [3 3]);
    verifySize(testCase, results.wpli_matrix, [3 3]);
    verifySize(testCase, results.plv_3d_matrix, [3 3 99]);
    verifySize(testCase, results.pli_3d_matrix, [3 3 99]);
    verifySize(testCase, results.wpli_3d_matrix, [3 3 99]);

    verifyEqual(testCase, diag(results.plv_matrix), ones(3, 1));
    verifyEqual(testCase, diag(results.pli_matrix), zeros(3, 1));
    verifyEqual(testCase, diag(results.wpli_matrix), zeros(3, 1));
    verifyEqual(testCase, results.plv_matrix, results.plv_matrix.', ...
        'AbsTol', 1e-14);
    verifyEqual(testCase, results.pli_matrix, results.pli_matrix.', ...
        'AbsTol', 1e-14);
    verifyEqual(testCase, results.wpli_matrix, results.wpli_matrix.', ...
        'AbsTol', 1e-14);

    verifyEqual(testCase, results.window_start_samples, 1:10:981);
    verifyEqual(testCase, results.window_length_samples, 20);
    verifyEqual(testCase, results.time_axis(1), 0.095, 'AbsTol', 1e-12);
    verifyEqual(testCase, results.time_axis(end), 9.895, 'AbsTol', 1e-12);

    all_values = [results.plv_matrix(:); results.pli_matrix(:); ...
        results.wpli_matrix(:); results.plv_3d_matrix(:); ...
        results.pli_3d_matrix(:); results.wpli_3d_matrix(:)];
    verifyTrue(testCase, all(isfinite(all_values)));
    verifyGreaterThanOrEqual(testCase, all_values, zeros(size(all_values)));
    verifyLessThanOrEqual(testCase, all_values, ones(size(all_values)));
end

function testRejectsInvalidInputs(testCase)
    unequal.ChannelA = zeros(100, 1);
    unequal.ChannelB = zeros(99, 1);
    verifyError(testCase, ...
        @() get_phase_connectivity(unequal, 100, 0.2, 0.1), ...
        'get_phase_connectivity:UnequalChannelLengths');

    short.ChannelA = zeros(10, 1);
    short.ChannelB = zeros(10, 1);
    verifyError(testCase, ...
        @() get_phase_connectivity(short, 100, 0.2, 0.1), ...
        'get_phase_connectivity:WindowExceedsSignal');
end
