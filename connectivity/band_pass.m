function filtered_eeg = band_pass(eeg_data, f_low, f_high, Fs, order)
%BAND_PASS Summary of this function goes here
%   Detailed explanation goes here

% 2. Design a Bandpass FIR filter
% We use 'bandpassfir' to create a filter object
bpFilter = designfilt('bandpassfir', ...
    'FilterOrder', order, ...      
    'CutoffFrequency1', f_low, ...
    'CutoffFrequency2', f_high, ...
    'SampleRate', Fs);

% 3. Apply filter to all 19 electrodes in the eeg_data structure
channel_names = fieldnames(eeg_data);
filtered_eeg = struct();

for i = 1:length(channel_names)
    ch_name = channel_names{i};
    raw_signal = eeg_data.(ch_name);
    
    % Use filtfilt for zero-phase distortion 
    filtered_eeg.(ch_name) = filtfilt(bpFilter, raw_signal);
end
end

