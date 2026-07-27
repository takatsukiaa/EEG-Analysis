clc;clearvars;close all;

filename = 'Raw_Data/No2.edf';

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

dummy_data.('Fp1_AVG') = 0;
dummy_data.('Fp2_AVG') = 0;
dummy_data.('F7_AVG') = 0;
dummy_data.('F3_AVG') = 0;
dummy_data.('FZ_AVG') = 0;
dummy_data.('F4_AVG') = 0;
dummy_data.('F8_AVG') = 0;
dummy_data.('T3_AVG') = 0;
dummy_data.('C3_AVG') = 0;
dummy_data.('CZ_AVG') = 0;
dummy_data.('C4_AVG') = 0;
dummy_data.('T4_AVG') = 0;
dummy_data.('T5_AVG') = 0;
dummy_data.('P3_AVG') = 0;
dummy_data.('PZ_AVG') = 0;
dummy_data.('P4_AVG') = 0;
dummy_data.('T6_AVG') = 0;
dummy_data.('O1_AVG') = 0;
dummy_data.('O2_AVG') = 0;

eeg_data = orderfields(eeg_data, dummy_data);
plot = 0;

%% No.1 Segment3 Pre
% StartTime3 = 118;
% EndTime3 = 123;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix);
% win_len = 0.5; % 1000 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment3 During
% StartTime3 = 123;
% EndTime3 = 124.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% 
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% --- 3. 計算兩兩電極間的滑動 PLV ---
% win_len = 0.5; % 500 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment3 Post 
% StartTime3 = 124.5;
% EndTime3 = 129.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix);
% % --- 3. 計算兩兩電極間的滑動 PLV ---
% win_len = 0.5; % 500 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment4 Pre 
% StartTime3 = 132;
% EndTime3 = 137;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment4 During 
% StartTime3 = 137;
% EndTime3 = 138;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% % --- 3. 計算兩兩電極間的滑動 PLV ---
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment4 Post 
% StartTime3 = 138;
% EndTime3 = 143;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix);
% --- 3. 計算兩兩電極間的滑動 PLV ---
% win_len = 0.2; % 500 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);


%% No.1 Segment5 Pre 
% StartTime3 = 139.5;
% EndTime3 = 144.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix);
% win_len = 0.5; % 500 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment5 During 
% StartTime3 = 144.5;
% EndTime3 = 146;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% % --- 3. 計算兩兩電極間的滑動 PLV ---
% win_len = 0.5; % 200 ms
% win_step = 0.1; % 100 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);

%% No.1 Segment5 Post 
% % StartTime3 = 146;
% % EndTime3 = 151;
% % 
% % eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% % 
% % filtered_eeg = band_pass(eeg_segment3, 25, 50, 200, 20);
% % 
% % phase_matrix = get_hilbert(filtered_eeg);
% % 
% % plv_matrix = get_section_PLV(phase_matrix);
% % --- 3. 計算兩兩電極間的滑動 PLV ---
% % win_len = 0.5; % 500 ms
% % win_step = 0.1; % 100 ms
% % fprintf('正在計算滑動視窗 PLV...\n');
% % [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step);


%% No.2 Segment1 Pre
% StartTime3 = 370.5;
% EndTime3 = 376.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% plot = 1;
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);

%% No.2 Segment1 During
% StartTime3 = 376.5;
% EndTime3 = 377.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% plot = 0;
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);


%% No.2 Segment1 Post
% StartTime3 = 383.5;
% EndTime3 = 389.5;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% plot = 1;
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);

%% No.2 Segment2 Pre
% StartTime3 = 1959;
% EndTime3 = 1964;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% plot = 1;
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);

%% No.2 Segment2 During
StartTime3 = 1964;
EndTime3 = 1965 ;

eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);

filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);

phase_matrix = get_hilbert(filtered_eeg);
plot = 0;

plv_matrix = get_section_PLV(phase_matrix, plot);
win_len = 0.2; % 200 ms
win_step = 0.1; % 500 ms
fprintf('正在計算滑動視窗 PLV...\n');
[plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);


%% No.2 Segment2 Post
% StartTime3 = 1969;
% EndTime3 = 1974;
% 
% eeg_segment3 = Segmentation(eeg_data, Fs, StartTime3, EndTime3);
% 
% filtered_eeg = band_pass(eeg_segment3, 35, 45, 200, 20);
% 
% phase_matrix = get_hilbert(filtered_eeg);
% plot = 1;
% 
% plv_matrix = get_section_PLV(phase_matrix, plot);
% win_len = 0.2; % 200 ms
% win_step = 0.1; % 500 ms
% fprintf('正在計算滑動視窗 PLV...\n');
% [plv_3d_matrix, t_axis] = get_sliding_PLV(filtered_eeg, Fs, win_len, win_step, plot);

