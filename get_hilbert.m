function phase_matrix = get_hilbert(filtered_eeg)
%HILBERT Summary of this function goes here
%   Detailed explanation goes here
channel_names = fieldnames(filtered_eeg);
for i = 1:length(channel_names)
    ch_name = channel_names{i};
    % 執行 Hilbert 變換並取得相位 (弧度)
    % angle 函數會回傳 -pi 到 pi 之間的瞬時相位
    phase_matrix.(ch_name) = angle(hilbert(filtered_eeg.(ch_name)));
end
end

