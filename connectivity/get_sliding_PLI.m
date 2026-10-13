function [pli_3d_matrix, wpli_3d_matrix, time_axis, channel_names] = ...
    get_sliding_PLI(filtered_eeg, Fs, win_sec, step_sec)
%GET_SLIDING_PLI Calculate sliding-window PLI and wPLI connectivity.
%   [PLI_3D, WPLI_3D, TIME_AXIS, CHANNEL_NAMES] = GET_SLIDING_PLI(...)
%   returns channel-by-channel-by-window matrices. TIME_AXIS contains each
%   window center in seconds relative to the beginning of FILTERED_EEG.
%
%   FILTERED_EEG must be a scalar structure of equally sized, real-valued
%   channel vectors and should already be band-pass filtered. The Hilbert
%   transform is calculated once for the complete segment before it is
%   divided into windows.

    validateattributes(Fs, {'numeric'}, {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'Fs', 2);
    validateattributes(win_sec, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'win_sec', 3);
    validateattributes(step_sec, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'step_sec', 4);

    [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg);
    [num_samples, num_channels] = size(analytic_signals);

    window_length = round(win_sec * Fs);
    window_step = round(step_sec * Fs);
    if window_length < 2
        error('get_sliding_PLI:WindowTooShort', ...
            'win_sec corresponds to fewer than two samples at Fs = %.6g Hz.', Fs);
    end
    if window_step < 1
        error('get_sliding_PLI:StepTooShort', ...
            'step_sec corresponds to fewer than one sample at Fs = %.6g Hz.', Fs);
    end
    if window_length > num_samples
        error('get_sliding_PLI:WindowExceedsSignal', ...
            'The window is %d samples, but the segment contains only %d samples.', ...
            window_length, num_samples);
    end

    start_indices = 1:window_step:(num_samples - window_length + 1);
    num_windows = numel(start_indices);
    time_axis = (start_indices - 1 + (window_length - 1) / 2) / Fs;

    pli_3d_matrix = zeros(num_channels, num_channels, num_windows);
    wpli_3d_matrix = zeros(num_channels, num_channels, num_windows);

    for window_index = 1:num_windows
        idx = start_indices(window_index):(start_indices(window_index) + window_length - 1);
        window_signals = analytic_signals(idx, :);

        for i = 1:num_channels
            for j = (i + 1):num_channels
                [pli_value, wpli_value] = pair_metrics( ...
                    window_signals(:, i), window_signals(:, j));

                pli_3d_matrix(i, j, window_index) = pli_value;
                pli_3d_matrix(j, i, window_index) = pli_value;
                wpli_3d_matrix(i, j, window_index) = wpli_value;
                wpli_3d_matrix(j, i, window_index) = wpli_value;
            end
        end
    end
end

function [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg)
    if ~isstruct(filtered_eeg) || ~isscalar(filtered_eeg)
        error('get_sliding_PLI:InvalidInput', ...
            'filtered_eeg must be a scalar structure containing channel vectors.');
    end

    channel_names = fieldnames(filtered_eeg);
    if numel(channel_names) < 2
        error('get_sliding_PLI:NotEnoughChannels', ...
            'At least two EEG channels are required to calculate connectivity.');
    end

    first_signal = filtered_eeg.(channel_names{1});
    validate_channel(first_signal, channel_names{1});
    num_samples = numel(first_signal);
    if num_samples < 2
        error('get_sliding_PLI:NotEnoughSamples', ...
            'Each EEG channel must contain at least two samples.');
    end

    signal_matrix = zeros(num_samples, numel(channel_names));
    for i = 1:numel(channel_names)
        signal = filtered_eeg.(channel_names{i});
        validate_channel(signal, channel_names{i});
        if numel(signal) ~= num_samples
            error('get_sliding_PLI:UnequalChannelLengths', ...
                'Channel "%s" has %d samples; expected %d.', ...
                channel_names{i}, numel(signal), num_samples);
        end
        signal_matrix(:, i) = double(signal(:));
    end

    analytic_signals = hilbert(signal_matrix);
end

function validate_channel(signal, channel_name)
    if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || isempty(signal)
        error('get_sliding_PLI:InvalidChannel', ...
            'Channel "%s" must be a nonempty, real numeric vector.', channel_name);
    end
    if any(~isfinite(signal(:)))
        error('get_sliding_PLI:NonfiniteChannel', ...
            'Channel "%s" contains NaN or Inf values.', channel_name);
    end
end

function [pli_value, wpli_value] = pair_metrics(z1, z2)
    imaginary_cross_spectrum = imag(z1 .* conj(z2));
    pli_value = abs(mean(sign(imaginary_cross_spectrum)));

    denominator = sum(abs(imaginary_cross_spectrum));
    if denominator == 0
        wpli_value = 0;
    else
        wpli_value = abs(sum(imaginary_cross_spectrum)) / denominator;
    end

    pli_value = min(1, max(0, pli_value));
    wpli_value = min(1, max(0, wpli_value));
end
