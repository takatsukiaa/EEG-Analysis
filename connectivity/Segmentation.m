function out_eeg = Segmentation(in_eeg, Fs, StartTime, EndTime)
%SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here
out_eeg = struct();
    channel_names = fieldnames(in_eeg);
    num_channels = length(channel_names);
    
    for i = 1:num_channels
        ch = channel_names{i};
        
        raw_vector = in_eeg.(ch);
        
        out_eeg.(ch) = raw_vector(Fs*StartTime:Fs*EndTime);
    end
end

