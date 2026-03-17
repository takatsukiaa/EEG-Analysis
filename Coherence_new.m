clc;clearvars;close all;

filename = 'Raw Data/No2.edf';

[~, ~, ext] = fileparts(filename);

%% File Reading
if strcmpi(ext, '.mat')
    disp(['偵測到 .mat 格式，正在載入：', filename]);
    eeg_data = load(filename);
    Fs = eeg_data.Fs;

    redundant_fields = {'data', 'Fs', 'C4_segment1', 'C4_segment2', 'O1_segment1', 'O1_segment2', 'O2_segment1', 'O2_segment2', 'P4_segment1', 'P4_segment2'};
    
    % fields_to_remove = intersect(fieldnames(eeg_data), redundant_fields);
    eeg_data = rmfield(eeg_data, redundant_fields);

elseif strcmpi(ext, '.edf')
    % --------------------------------------
    % 處理 .edf 檔案
    % --------------------------------------
    disp(['偵測到 .edf 格式，正在載入：', filename]);
    info = edfinfo(filename);
    
    fs = info.NumSamples ./ seconds(info.DataRecordDuration);
    Fs = fs(1,1);
    
    edf_table = edfread(filename);
    channel_names = edf_table.Properties.VariableNames;
    
    for i = 1:length(channel_names)
        ch_name = channel_names{i};
        
        % 忽略非腦波的通道
        if contains(ch_name, 'emg', 'IgnoreCase', true) || contains(ch_name, 'ekg', 'IgnoreCase', true)
            continue;
        end
        
        % 將分塊的數據垂直串接 (vertcat) 攤平
        eeg_data.(ch_name) = vertcat(edf_table.(ch_name){:});
    end

else
    error(['不支援的檔案格式：', ext, '！請提供 .mat 或 .edf 檔案。']);
end

disp(['資料載入完成！採樣率 Fs = ', num2str(Fs), ' Hz，共 ', num2str(length(fieldnames(eeg_data))), ' 個通道。']);

%% Segment1
StartTime1 = 1960;
EndTime1 = 1970;

eeg_segment1 = Segmentation(eeg_data, Fs, StartTime1, EndTime1);

output_file1 = 'Results/No.2 coherence_segment1.txt';
analyze_envelope_xcorr(eeg_segment1, Fs, output_file1);

%% Segment2
StartTime2 = 370;
EndTime2 = 390;

eeg_segment2 = Segmentation(eeg_data, Fs, StartTime2, EndTime2);

output_file2 = 'Results/No.2 coherence_segment2.txt';
analyze_envelope_xcorr(eeg_segment2, Fs, output_file2);

%% Segment3
% StartTime3 = 122;
% EndTime3 = 126;
% 
% eeg_segment2 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% output_file2 = 'Results/No.1 coherence_segment3.txt';
% analyze_envelope_xcorr(eeg_segment2, Fs, output_file2);

%% Segment4
% StartTime3 = 135;
% EndTime3 = 140;
% 
% eeg_segment2 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% output_file2 = 'Results/No.1 coherence_segment4.txt';
% analyze_envelope_xcorr(eeg_segment2, Fs, output_file2);

%% Segment5
% StartTime3 = 143;
% EndTime3 = 147;
% 
% eeg_segment2 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% output_file2 = 'Results/No.1 coherence_segment5.txt';
% analyze_envelope_xcorr(eeg_segment2, Fs, output_file2);



