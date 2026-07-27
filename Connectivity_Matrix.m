clc; clearvars; close all;

% 1. 讀取剛才生成的 txt 檔案
% 請確保檔名與您產生的檔案名稱一致
filename = 'Results/No.1 coherence_segment3_post.txt';
T = readtable(filename, 'Delimiter', '\t');

% 2. 抓取所有不重複的通道名稱，以建立矩陣的維度
% 將 Channel_A 和 Channel_B 合併後取唯一值
all_channels = unique([T.Channel_A; T.Channel_B], 'stable');
num_ch = length(all_channels);

% 3. 初始化一個 N x N 的矩陣
% 對角線 (自己對自己) 的相關係數預設為 1
R_matrix = eye(num_ch); 
Lag_matrix = zeros(num_ch); 

% 4. 將 txt 檔內的數據填入矩陣中
for i = 1:height(T)
    chA = T.Channel_A{i};
    chB = T.Channel_B{i};
    r_val = T.R_value(i);
    lag_val = T.Lag_point(i);
    
    % 找出這兩個通道在矩陣中的索引位置
    idxA = find(strcmp(all_channels, chA));
    idxB = find(strcmp(all_channels, chB));
    
    % 填入 R 值 (無方向性，所以矩陣是對稱的)
    R_matrix(idxA, idxB) = r_val;
    R_matrix(idxB, idxA) = r_val; 
    
    % 填入 Lag 值 (有方向性，所以矩陣是反對稱的，A相對於B為正，B相對於A就為負)
    Lag_matrix(idxA, idxB) = lag_val;
    Lag_matrix(idxB, idxA) = -lag_val; 
end

% ==========================================
% 5. 繪製 R 值 (相關性) 的熱力圖
% ==========================================
figure('Name', 'Connectivity Matrix (R-value)', 'Position', [100, 100, 800, 650]);

% 使用 MATLAB 內建的 heatmap 函數
h1 = heatmap(all_channels, all_channels, R_matrix);

% 外觀設定
h1.Title = 'Slow-Gamma Connectivity Matrix (R-value)';
h1.XLabel = 'Electrodes';
h1.YLabel = 'Electrodes';

% 設定顏色範圍 (Colormap)
% 為了讓對比更明顯，我們將最低顏色設為數據中的最小值，最高設為 1
min_R = min(T.R_value);
h1.ColorLimits = [min_R, 1]; 

% 更換為大腦連結常用的 'jet' 或 'parula' 色系
colormap('jet'); 

% ==========================================
% 6. (選用) 繪製 Lag 值 (延遲) 的熱力圖
% ==========================================
figure('Name', 'Lag Matrix (Points)', 'Position', [950, 100, 800, 650]);

h2 = heatmap(all_channels, all_channels, Lag_matrix);
h2.Title = 'Slow-Gamma Lag Matrix (Points)';
h2.XLabel = 'Electrodes';
h2.YLabel = 'Electrodes';

% Lag 通常有正有負，我們將顏色對稱於 0 設定
max_lag = max(abs(T.Lag_point));
if max_lag == 0
    max_lag = 1; % 避免除以零
end
h2.ColorLimits = [-max_lag, max_lag];
colormap(h2, 'hot'); % 可以用不同色系區隔