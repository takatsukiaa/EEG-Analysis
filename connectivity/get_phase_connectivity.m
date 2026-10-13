function results = get_phase_connectivity(filtered_eeg, Fs, win_sec, step_sec)
%GET_PHASE_CONNECTIVITY Calculate aligned PLV, PLI, and wPLI connectivity.
%   RESULTS = GET_PHASE_CONNECTIVITY(FILTERED_EEG, FS, WIN_SEC, STEP_SEC)
%   calculates whole-segment and sliding-window connectivity from the same
%   analytic signals and sample windows. FILTERED_EEG must be a scalar
%   structure of equally sized, band-pass-filtered channel vectors.
%
%   RESULTS contains:
%       plv_matrix, pli_matrix, wpli_matrix
%       plv_3d_matrix, pli_3d_matrix, wpli_3d_matrix
%       time_axis, channel_names, window_start_samples, window_length_samples
%
%   PLV has a diagonal of one. PLI and wPLI have diagonals of zero because
%   self-connectivity is not a lagged connection.

    validateattributes(Fs, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'Fs', 2);
    validateattributes(win_sec, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'win_sec', 3);
    validateattributes(step_sec, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'step_sec', 4);

    [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg);
    [num_samples, num_channels] = size(analytic_signals);

    window_length = round(win_sec * Fs);
    window_step = round(step_sec * Fs);
    if window_length < 2
        error('get_phase_connectivity:WindowTooShort', ...
            'win_sec corresponds to fewer than two samples at Fs = %.6g Hz.', Fs);
    end
    if window_step < 1
        error('get_phase_connectivity:StepTooShort', ...
            'step_sec corresponds to fewer than one sample at Fs = %.6g Hz.', Fs);
    end
    if window_length > num_samples
        error('get_phase_connectivity:WindowExceedsSignal', ...
            'The window is %d samples, but the segment contains only %d samples.', ...
            window_length, num_samples);
    end

    start_indices = 1:window_step:(num_samples - window_length + 1);
    num_windows = numel(start_indices);
    time_axis = (start_indices - 1 + (window_length - 1) / 2) / Fs;

    [plv_matrix, pli_matrix, wpli_matrix] = ...
        calculate_matrices(analytic_signals);

    plv_3d_matrix = zeros(num_channels, num_channels, num_windows);
    pli_3d_matrix = zeros(num_channels, num_channels, num_windows);
    wpli_3d_matrix = zeros(num_channels, num_channels, num_windows);

    for window_index = 1:num_windows
        idx = start_indices(window_index):(start_indices(window_index) + window_length - 1);
        [plv_3d_matrix(:, :, window_index), ...
         pli_3d_matrix(:, :, window_index), ...
         wpli_3d_matrix(:, :, window_index)] = ...
            calculate_matrices(analytic_signals(idx, :));
    end

    results = struct( ...
        'plv_matrix', plv_matrix, ...
        'pli_matrix', pli_matrix, ...
        'wpli_matrix', wpli_matrix, ...
        'plv_3d_matrix', plv_3d_matrix, ...
        'pli_3d_matrix', pli_3d_matrix, ...
        'wpli_3d_matrix', wpli_3d_matrix, ...
        'time_axis', time_axis, ...
        'channel_names', {channel_names}, ...
        'window_start_samples', start_indices, ...
        'window_length_samples', window_length);
end

function [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg)
    if ~isstruct(filtered_eeg) || ~isscalar(filtered_eeg)
        error('get_phase_connectivity:InvalidInput', ...
            'filtered_eeg must be a scalar structure containing channel vectors.');
    end

    channel_names = fieldnames(filtered_eeg);
    if numel(channel_names) < 2
        error('get_phase_connectivity:NotEnoughChannels', ...
            'At least two EEG channels are required to calculate connectivity.');
    end

    first_signal = filtered_eeg.(channel_names{1});
    validate_channel(first_signal, channel_names{1});
    num_samples = numel(first_signal);
    if num_samples < 2
        error('get_phase_connectivity:NotEnoughSamples', ...
            'Each EEG channel must contain at least two samples.');
    end

    signal_matrix = zeros(num_samples, numel(channel_names));
    for channel_index = 1:numel(channel_names)
        channel_name = channel_names{channel_index};
        signal = filtered_eeg.(channel_name);
        validate_channel(signal, channel_name);
        if numel(signal) ~= num_samples
            error('get_phase_connectivity:UnequalChannelLengths', ...
                'Channel "%s" has %d samples; expected %d.', ...
                channel_name, numel(signal), num_samples);
        end
        signal_matrix(:, channel_index) = double(signal(:));
    end

    % The one Hilbert transform for this analysis is performed here, before
    % any windowing, so all metrics share identical analytic samples and do
    % not acquire separate short-window edge artifacts.
    analytic_signals = hilbert(signal_matrix);
end

function validate_channel(signal, channel_name)
    if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || isempty(signal)
        error('get_phase_connectivity:InvalidChannel', ...
            'Channel "%s" must be a nonempty, real numeric vector.', channel_name);
    end
    if any(~isfinite(signal(:)))
        error('get_phase_connectivity:NonfiniteChannel', ...
            'Channel "%s" contains NaN or Inf values.', channel_name);
    end
end

function [plv_matrix, pli_matrix, wpli_matrix] = calculate_matrices(analytic_signals)
    num_channels = size(analytic_signals, 2);
    plv_matrix = eye(num_channels);
    pli_matrix = zeros(num_channels, num_channels);
    wpli_matrix = zeros(num_channels, num_channels);

    for i = 1:num_channels
        for j = (i + 1):num_channels
            [plv_value, pli_value, wpli_value] = pair_metrics( ...
                analytic_signals(:, i), analytic_signals(:, j));

            plv_matrix(i, j) = plv_value;
            plv_matrix(j, i) = plv_value;
            pli_matrix(i, j) = pli_value;
            pli_matrix(j, i) = pli_value;
            wpli_matrix(i, j) = wpli_value;
            wpli_matrix(j, i) = wpli_value;
        end
    end
end

function [plv_value, pli_value, wpli_value] = pair_metrics(z1, z2)
    normalized_phase_difference = (z1 ./ abs(z1)) .* conj(z2 ./ abs(z2));
    invalid_phase = ~isfinite(normalized_phase_difference);
    normalized_phase_difference(invalid_phase) = 0;
    valid_count = numel(normalized_phase_difference) - nnz(invalid_phase);
    if valid_count == 0
        plv_value = 0;
    else
        plv_value = abs(sum(normalized_phase_difference)) / valid_count;
    end

    imaginary_cross_spectrum = imag(z1 .* conj(z2));
    pli_value = abs(mean(sign(imaginary_cross_spectrum)));

    denominator = sum(abs(imaginary_cross_spectrum));
    if denominator == 0
        wpli_value = 0;
    else
        wpli_value = abs(sum(imaginary_cross_spectrum)) / denominator;
    end

    plv_value = min(1, max(0, plv_value));
    pli_value = min(1, max(0, pli_value));
    wpli_value = min(1, max(0, wpli_value));
end
