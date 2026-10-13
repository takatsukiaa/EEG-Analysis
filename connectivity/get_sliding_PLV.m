function [plv_3d, time_axis] = get_sliding_PLV(data_struct, fs, win_sec, step_sec, plot)
% GET_SLIDING_PLV Calculates the sliding-window Phase Locking Value (PLV) 
% for multi-channel EEG data stored in a structure.
%
% Inputs:
%   data_struct - Structure containing filtered EEG signals (each field is a channel)
%   fs          - Sampling rate (Hz)
%   win_sec     - Window length in seconds
%   step_sec    - Window sliding step in seconds
%
% Outputs:
%   plv_3d      - 3D matrix (channels x channels x windows)
%   time_axis   - Vector of time points for each window center (seconds)

    % 1. Extract channel information to determine matrix size
    chan_names = fieldnames(data_struct);
    num_ch = length(chan_names);
    num_samples = length(data_struct.(chan_names{1}));
    
    % 2. Extract instantaneous phases using Hilbert Transform
    phase_matrix = zeros(num_ch, num_samples);
    for i = 1:num_ch
        phase_matrix(i, :) = angle(hilbert(data_struct.(chan_names{i})));
    end

    % 3. Setup window parameters in samples
    win_len = round(win_sec * fs);
    win_step = round(step_sec * fs);
    start_indices = 1:win_step:(num_samples - win_len + 1);
    num_windows = length(start_indices);
    
    % 4. Pre-allocate results and compute time axis (center of each window)
    plv_3d = zeros(num_ch, num_ch, num_windows);
    time_axis = (start_indices + (win_len / 2) - 1) / fs;

    % 5. Compute PLV for each window pair-wise
    for w = 1:num_windows
        idx = start_indices(w) : (start_indices(w) + win_len - 1);
        win_phases = phase_matrix(:, idx);
        
        for i = 1:num_ch
            for j = i:num_ch
                if i == j
                    plv_3d(i, j, w) = 1;
                else
                    % Compute phase difference complex vector
                    phase_diff = win_phases(i, :) - win_phases(j, :);
                    % PLV is the amplitude of the mean phase difference vector
                    val = abs(sum(exp(1i * phase_diff))) / win_len;
                    
                    plv_3d(i, j, w) = val;
                    plv_3d(j, i, w) = val; % Symmetric property
                end
            end
        end
    end
    fprintf('計算完成！\n');
    


    if plot
        output_folder = 'PLV_Results_Images/No2_segment2_post_slice_5';
        if ~exist(output_folder, 'dir'), mkdir(output_folder); end
        
        num_windows = size(plv_3d, 3);
        
        for w = 1:num_windows
        % 1. 建立隱藏視窗並設定大尺寸
            fig = figure('Visible', 'off', 'Color', 'w');
            set(fig, 'Position', [0, 0, 1200, 1200]); 

            current_plv = plv_3d(:,:,w);
            mask = triu(true(size(current_plv)), 1); % Logical mask for strictly upper triangle
            current_plv(mask) = NaN;
        
            % 2. 繪製熱圖
            h = heatmap(chan_names, chan_names, current_plv);
            h.Title = sprintf('Window %d - Time: %.2f s', w, time_axis(w));
            h.Colormap = jet;
            h.ColorLimits = [0 1];
            h.CellLabelFormat = '%.2f';
            h.FontSize = 8; % 調整字體大小
            h.MissingDataColor = [1 1 1];
        
            % 3. 高解析度存檔 (300 DPI)
            filename = fullfile(output_folder, sprintf('PLV_Slice5_Window_%03d.png', w));
            exportgraphics(fig, filename, 'Resolution', 300);
        
            close(fig); % 關閉以釋放記憶體
        end
        fprintf('所有視窗已儲存至 %s 資料夾。\n', output_folder);
    end
end