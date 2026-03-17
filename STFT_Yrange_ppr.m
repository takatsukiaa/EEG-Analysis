clear; clc;
close all;

%% 1. 讀取 EDF 資料
[file, path] = uigetfile('*.edf', '選擇腦波 EDF 檔案');
if isequal(file, 0), return; end
fullPath = fullfile(path, file);

edfData = edfread(fullPath);
data_full = cell2mat(table2cell(edfData))'; % 先讀取完整資料

% 取得採樣率
info = edfinfo(fullPath);
fs = info.NumSamples(1) / seconds(info.DataRecordDuration); 

%% 2. --- 新增：截取特定 60 秒 ---
startTime = 954;  % <--- 在此設定起點（秒），例如從第 30 秒開始
duration  = 90;  % 設定持續時間為 60 秒

startIdx = round(startTime * fs) + 1;
endIdx   = startIdx + round(duration * fs) - 1;

% 檢查索引是否超出檔案範圍
if endIdx > size(data_full, 2)
    fprintf('警告：要求的時間範圍超出檔案長度，將擷取至檔案末尾。\n');
    endIdx = size(data_full, 2);
end

% 擷取資料
data = data_full(:, startIdx:endIdx);
data(20:21, :) = [];

% 重新定義時間軸（從 0 開始計時，或保留原始時間起點）
% 如果要圖表顯示 0~60 秒：
t = (0:size(data,2)-1) / fs;
% 原有的通道名稱與時間軸定義
channelNames = edfData.Properties.VariableNames;
channelNames(:,20:21) = [];

%% 參數設定 (調整以觀察 Gamma)
innerMargin = 0.025;      
intra_gap   = 0.005;      

windowLen = 2*fs;
overlapLen = fs;
nfft = 2^nextpow2(windowLen);
cmap = jet(128);  

% --- 修改點：擴展頻段至 Gamma (30-80Hz) ---
low_cutoff = 30;
high_cutoff = 99; % 調整為 80Hz 以涵蓋 Gamma
filter_order = 4;
Wn = [low_cutoff, high_cutoff] / (fs/2);
[b, a] = butter(filter_order, Wn, 'bandpass');

% --- 新增點：60Hz 陷波濾波器 (濾除電力線雜訊) ---
[b_n, a_n] = iirnotch(60/(fs/2), (60/(fs/2))/35);

% 通道分組與圖表標題 (維持原樣)
channel_groups = {3:6, 1:2, 11:14, 15:16, 7:10, [1, 2, 5, 6, 9, 10], [1, 2, 5, 6], [5, 6, 9, 10], [3, 4, 5, 6, 9, 10]};
titles = {'EEG Signals & Spectrogram - F Channels', ...
    'EEG Signals & Spectrogram - Fp Channels', ...
    'EEG Signals & Spectrogram - C & P Channels', ...
    'EEG Signals & Spectrogram - O Channels', ...
    'EEG Signals & Spectrogram - T Channels', ...
    'EEG Signals & Spectrogram - Fp Channels & C Channels, O Channels', ...
    'EEG Signals & Spectrogram - Fp Channels & C Channels, O Channels', ...
    'EEG Signals & Spectrogram - Fp Channels & C Channels, O Channels', ...
    'EEG Signals & Spectrogram - F Channels & C Channels, O Channels'};
layout_settings = {[2,2], [2,2], [2,2], [2, 2], [2, 2], [2, 2], [2, 2], [2, 2],[3, 2]};

%% ======================
%% 時頻圖：百分位數顯示 (調整 Y 軸以觀察 Gamma)
%% ======================
onset_time = 10; 
offset_time = 45;
for figIdx = 1:5
    chIndices = channel_groups{figIdx};
    nRows = layout_settings{figIdx}(1);
    nCols = layout_settings{figIdx}(2);
    
    all_power_values = [];
    for i = 1:length(chIndices)
        ch = chIndices(i);
        % 套用帶通與陷波濾波
        eegSignal = filtfilt(b, a, data(ch, :));
        eegSignal = filtfilt(b_n, a_n, eegSignal); 
        [s, ~, ~] = spectrogram(eegSignal, windowLen, overlapLen, nfft, fs);
        sPow = abs(s).^2;       
        all_power_values = [all_power_values; sPow(:)];
    end
    
    figure('Position',[100, 50, 1600, 900]);
    sgtitle([titles{figIdx} ' (Gamma Range)'], 'FontSize', 16);
    
    % 佈局參數設定 (與原碼一致)
    leftMargin = 0.1; rightMargin = 0.1;
    topMargin = 0.1; bottomMargin = 0.1;
    spacingX = 0.03; spacingY = 0.03;
    tileWidth = (1 - leftMargin - rightMargin - (nCols-1)*spacingX) / nCols;
    tileHeight = (1 - topMargin - bottomMargin - (nRows-1)*spacingY) / nRows;
    
    specAxes = [];
    
    for i = 1:length(chIndices)
        ch = chIndices(i);
        row = ceil(i / nCols);
        col = mod(i-1, nCols) + 1;
        posX = leftMargin + (col-1)*(tileWidth+spacingX);
        posY = bottomMargin + (nRows - row)*(tileHeight+spacingY);
        
        topHeight = (tileHeight - intra_gap)/2;
        axTime_pos = [posX + innerMargin, posY + topHeight + 0.5*intra_gap + innerMargin, tileWidth - 2*innerMargin, topHeight - innerMargin];
        axSpec_pos = [posX + innerMargin, posY + innerMargin, tileWidth - 2*innerMargin, topHeight - innerMargin];
        
        % 時域圖
        axTime = axes('Position', axTime_pos);
        eegSignal = filtfilt(b, a, data(ch, :));
        eegSignal = filtfilt(b_n, a_n, eegSignal);
        plot(axTime, t, eegSignal, 'b');
        title(axTime, channelNames{ch}, 'Interpreter', 'none', 'FontSize', 10);
        xlim(axTime, [t(1) t(end)]);
        ylim(axTime, [-100 100]);
        set(axTime, 'XTick', []);
        
        % 頻譜圖 (調整 Y 軸)
        axSpec = axes('Position', axSpec_pos);
        specAxes = [specAxes, axSpec];
        [s, f, t_s] = spectrogram(eegSignal, windowLen, overlapLen, nfft, fs);
        sPow = abs(s).^2;
        
        freq_range = (f >= low_cutoff & f <= high_cutoff);
        sPlot_linear = sPow(freq_range, :);
        fPlot = f(freq_range);
        
        lowerBound = prctile(all_power_values, 0);
        upperBound = prctile(all_power_values, 99);
        sPlot_percentile = (sPlot_linear - lowerBound) / (upperBound - lowerBound) * 100;
        sPlot_percentile(sPlot_percentile < 0) = 0;
        sPlot_percentile(sPlot_percentile > 100) = 100;
        
        imagesc(axSpec, t_s, fPlot, sPlot_percentile);
        axis(axSpec, 'xy');
        caxis(axSpec, [0, 100]);         
        colormap(axSpec, cmap);
        hold(axSpec, 'on'); % 確保標記線畫在熱圖上方

        % 畫垂直線：'r--' 代表紅色虛線
        xl_start = xline(axSpec, onset_time, 'r--', 'LineWidth', 2);
        xl_end   = xline(axSpec, offset_time, 'r--', 'LineWidth', 2);
        
        % 加入文字標註 (放在標記線上方)
        text(axSpec, onset_time, 75, ' Onset', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);
        text(axSpec, offset_time, 75, ' Offset', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);


        % --- 修改點：Y 軸顯示到 80Hz ---
        ylim(axSpec, [30 80]);
        yticks(axSpec, [30, 40, 50, 60, 70, 80]);
        xlabel(axSpec, 'Time (SEC)', 'FontSize', 8);
        ylabel(axSpec, 'Freq (Hz)', 'FontSize', 8);
    end
    
    if ~isempty(specAxes)
        cb = colorbar(specAxes(end));
        posSpec = get(specAxes(end), 'Position');
        cb.Position = [posSpec(1) + posSpec(3) + 0.01, posSpec(2), 0.02, posSpec(4)];
        ylabel(cb, 'Power %', 'FontSize', 8);
    end
end

%% ================== PSD 繪圖 (調整 X 軸) ==================
figure('Name','PSD Gamma','Position',[200 50 1000 800]);
% 重新定義分組 (根據您的 sChannel 順序)
% F Channels: 3, 4, 5, 6
% O Channels: 11, 12
% C Channels: 15, 16
% T Channels: 7, 8, 9, 10
group_indices = { [3, 4, 5, 6], [1, 2], [11, 12], [15, 16], [7,8,9,10], [13,14] }; 
group_names   = {"F Channels", "Fp Channels", "C Channels", "O Channels", "T Channels","P Channels"};
% 對應的 Legend
group_legends = { ["F4_AVG", "F3_AVG", "F8_AVG", "F7_AVG"], ...
                  ["Fp2_AVG","Fp1_AVG"], ...
                  ["C4_AVG", "C3_AVG"], ...
                  ["O2_AVG", "O1_AVG"] , ...
                  ["T4_AVG", "T3_AVG", "T6_AVG", "T5_AVG"], ...
                  ["P4_AVG", "P3_AVG"]};

% 設定顏色順序 (確保同一圖內不同曲線好辨識)
colors = ['b', 'r', 'g', 'm', 'c', 'k']; 

for ch = 1:length(group_indices)
    subplot(3, 2, ch);
    hold on; grid on;
    
    current_channels = group_indices{ch};
    current_legends  = group_legends{ch};
    
    % 在同一個子圖內遍歷該組的所有通道
    for i = 1:length(current_channels)
        target_ch = current_channels(i);
        
        % 訊號處理 (帶通 + 陷波)
        sig = filtfilt(b_n, a_n, filtfilt(b, a, data(target_ch, :)));
        
        % 計算 PSD
        [pxx, fxx] = pwelch(sig, windowLen, overlapLen, nfft, fs);
        
        % 繪圖，循環使用顏色
        plot(fxx, pxx, 'Color', colors(mod(i-1,6)+1), 'LineWidth', 1);
    end
    
    % 圖表美化
    legend(current_legends);
    xlim([0.5 99]);          % 聚焦 Gamma Band
    ylim([0 4]);           % 調整 Y 軸以觀察微弱的 Gamma 訊號
    xlabel('Frequency (Hz)');
    ylabel('Power (\muV²/Hz)');
    title(group_names{ch});
end