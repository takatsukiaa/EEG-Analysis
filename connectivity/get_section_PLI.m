function [pli_matrix, wpli_matrix, channel_names] = get_section_PLI(filtered_eeg)
%GET_SECTION_PLI Calculate whole-segment PLI and wPLI connectivity matrices.
%   [PLI_MATRIX, WPLI_MATRIX, CHANNEL_NAMES] = GET_SECTION_PLI(FILTERED_EEG)
%   accepts a scalar structure whose fields are equally sized, real-valued
%   channel vectors. FILTERED_EEG should already be band-pass filtered.
%
%   PLI follows Stam et al. (2007):
%       abs(mean(sign(imag(z1 .* conj(z2)))))
%   wPLI is the amplitude-weighted version:
%       abs(sum(imagCross)) / sum(abs(imagCross))
%   where z1 and z2 are Hilbert analytic signals.
%
%   Both output matrices are symmetric, bounded by [0, 1], and have zero
%   diagonals because self-connectivity is not a lagged connection.

    [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg);
    num_channels = size(analytic_signals, 2);

    pli_matrix = zeros(num_channels, num_channels);
    wpli_matrix = zeros(num_channels, num_channels);

    for i = 1:num_channels
        for j = (i + 1):num_channels
            [pli_value, wpli_value] = pair_metrics( ...
                analytic_signals(:, i), analytic_signals(:, j));

            pli_matrix(i, j) = pli_value;
            pli_matrix(j, i) = pli_value;
            wpli_matrix(i, j) = wpli_value;
            wpli_matrix(j, i) = wpli_value;
        end
    end
end

function [analytic_signals, channel_names] = prepare_analytic_signals(filtered_eeg)
    if ~isstruct(filtered_eeg) || ~isscalar(filtered_eeg)
        error('get_section_PLI:InvalidInput', ...
            'filtered_eeg must be a scalar structure containing channel vectors.');
    end

    channel_names = fieldnames(filtered_eeg);
    if numel(channel_names) < 2
        error('get_section_PLI:NotEnoughChannels', ...
            'At least two EEG channels are required to calculate connectivity.');
    end

    first_signal = filtered_eeg.(channel_names{1});
    validate_channel(first_signal, channel_names{1});
    num_samples = numel(first_signal);
    if num_samples < 2
        error('get_section_PLI:NotEnoughSamples', ...
            'Each EEG channel must contain at least two samples.');
    end

    signal_matrix = zeros(num_samples, numel(channel_names));
    for i = 1:numel(channel_names)
        signal = filtered_eeg.(channel_names{i});
        validate_channel(signal, channel_names{i});
        if numel(signal) ~= num_samples
            error('get_section_PLI:UnequalChannelLengths', ...
                'Channel "%s" has %d samples; expected %d.', ...
                channel_names{i}, numel(signal), num_samples);
        end
        signal_matrix(:, i) = double(signal(:));
    end

    % Hilbert is applied once to the complete segment. Sliding-window code
    % subsequently slices these analytic signals, avoiding window-edge
    % artifacts from repeatedly applying Hilbert to short windows.
    analytic_signals = hilbert(signal_matrix);
end

function validate_channel(signal, channel_name)
    if ~isnumeric(signal) || ~isreal(signal) || ~isvector(signal) || isempty(signal)
        error('get_section_PLI:InvalidChannel', ...
            'Channel "%s" must be a nonempty, real numeric vector.', channel_name);
    end
    if any(~isfinite(signal(:)))
        error('get_section_PLI:NonfiniteChannel', ...
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

    % Protect the public contract from negligible floating-point overshoot.
    pli_value = min(1, max(0, pli_value));
    wpli_value = min(1, max(0, wpli_value));
end
