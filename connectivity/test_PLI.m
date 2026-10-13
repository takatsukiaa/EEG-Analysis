function tests = test_PLI
%TEST_PLI Synthetic tests for get_section_PLI and get_sliding_PLI.
%   Run with:
%       results = runtests('test_PLI.m');
%       assertSuccess(results);

    tests = functiontests(localfunctions);
end

function testConstantNonzeroPhaseLag(testCase)
    Fs = 500;
    time = (0:(10 * Fs - 1)).' / Fs;
    signals.ChannelA = cos(2 * pi * 12 * time);
    signals.ChannelB = cos(2 * pi * 12 * time + pi / 4);

    [pli_matrix, wpli_matrix] = get_section_PLI(signals);

    verifyGreaterThan(testCase, pli_matrix(1, 2), 0.99);
    verifyGreaterThan(testCase, wpli_matrix(1, 2), 0.99);
end

function testIdenticalZeroLagSignals(testCase)
    Fs = 500;
    time = (0:(4 * Fs - 1)).' / Fs;
    signal = cos(2 * pi * 10 * time) + 0.2 * cos(2 * pi * 18 * time);
    signals.ChannelA = signal;
    signals.ChannelB = signal;

    [pli_matrix, wpli_matrix] = get_section_PLI(signals);

    verifyEqual(testCase, pli_matrix(1, 2), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, wpli_matrix(1, 2), 0, 'AbsTol', 1e-14);
end

function testIndependentSignalsAreNearZero(testCase)
    rng(17, 'twister');
    num_samples = 50000;
    signals.ChannelA = randn(num_samples, 1);
    signals.ChannelB = randn(num_samples, 1);

    [pli_matrix, wpli_matrix] = get_section_PLI(signals);

    verifyLessThan(testCase, pli_matrix(1, 2), 0.08);
    verifyLessThan(testCase, wpli_matrix(1, 2), 0.08);
end

function testMatrixContract(testCase)
    Fs = 200;
    time = (0:(5 * Fs - 1)).' / Fs;
    signals.C3 = cos(2 * pi * 8 * time);
    signals.C4 = cos(2 * pi * 8 * time + pi / 3);
    signals.Pz = cos(2 * pi * 8 * time - pi / 5);

    [pli_matrix, wpli_matrix, channel_names] = get_section_PLI(signals);

    verifyEqual(testCase, channel_names, {'C3'; 'C4'; 'Pz'});
    verifyEqual(testCase, pli_matrix, pli_matrix.', 'AbsTol', 1e-14);
    verifyEqual(testCase, wpli_matrix, wpli_matrix.', 'AbsTol', 1e-14);
    verifyEqual(testCase, diag(pli_matrix), zeros(3, 1));
    verifyEqual(testCase, diag(wpli_matrix), zeros(3, 1));
    verifyTrue(testCase, all(isfinite(pli_matrix), 'all'));
    verifyTrue(testCase, all(isfinite(wpli_matrix), 'all'));
    verifyGreaterThanOrEqual(testCase, pli_matrix, zeros(3));
    verifyLessThanOrEqual(testCase, pli_matrix, ones(3));
    verifyGreaterThanOrEqual(testCase, wpli_matrix, zeros(3));
    verifyLessThanOrEqual(testCase, wpli_matrix, ones(3));
end

function testSlidingDimensionsAndTimeAxis(testCase)
    Fs = 100;
    num_samples = 1000;
    time = (0:(num_samples - 1)).' / Fs;
    signals.ChannelA = cos(2 * pi * 7 * time);
    signals.ChannelB = cos(2 * pi * 7 * time + pi / 4);
    signals.ChannelC = cos(2 * pi * 9 * time);

    [pli_3d, wpli_3d, time_axis, channel_names] = ...
        get_sliding_PLI(signals, Fs, 0.2, 0.1);

    verifySize(testCase, pli_3d, [3 3 99]);
    verifySize(testCase, wpli_3d, [3 3 99]);
    verifyEqual(testCase, channel_names, ...
        {'ChannelA'; 'ChannelB'; 'ChannelC'});
    verifyEqual(testCase, time_axis(1), 0.095, 'AbsTol', 1e-12);
    verifyEqual(testCase, time_axis(end), 9.895, 'AbsTol', 1e-12);
    verifyEqual(testCase, squeeze(pli_3d(1, 2, :)), ...
        squeeze(pli_3d(2, 1, :)), 'AbsTol', 1e-14);
    verifyEqual(testCase, squeeze(wpli_3d(1, 2, :)), ...
        squeeze(wpli_3d(2, 1, :)), 'AbsTol', 1e-14);
end

function testRejectsUnequalChannelLengths(testCase)
    signals.ChannelA = zeros(100, 1);
    signals.ChannelB = zeros(99, 1);

    verifyError(testCase, @() get_section_PLI(signals), ...
        'get_section_PLI:UnequalChannelLengths');
    verifyError(testCase, @() get_sliding_PLI(signals, 100, 0.2, 0.1), ...
        'get_sliding_PLI:UnequalChannelLengths');
end
