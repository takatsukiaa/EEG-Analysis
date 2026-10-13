clear; clc;
close all;

fileName = 'No.1 Sz.edf';
[data, anns] = edfread(fileName); 
info = edfinfo(fileName)

% --- 1. 參數設定 ---
targetChs1 = {'Fp1_AVG', 'Fp2_AVG', 'F3_AVG', 'F4_AVG', 'F7_AVG', 'F8_AVG', 'T3_AVG', 'T4_AVG',};
targetChs2 = {'T5_AVG', 'T6_AVG', 'C3_AVG', 'C4_AVG', 'P3_AVG', 'P4_AVG', 'O1_AVG', 'O2_AVG'};

fs = 200;                   % 採樣率
secPerRow = 7;              % 每個 Row 的秒數
startTime = 954;            % 起始秒數
duration = 90;              % 擷取長度 (秒)

% --- 2. 索引計算 ---
startRow = floor(startTime / secPerRow) + 1;
endRow = ceil((startTime + duration) / secPerRow);
numSamples = duration * fs;

% --- 3. 批次擷取與繪圖 ---
figure('Color', 'w', 'Name', 'Raw EEG Signal (954s-1014s)');
timeVec = linspace(startTime, startTime + duration, numSamples);

yLimit = [-200, 200];
ax = zeros(1, length(targetChs1)); % 儲存每個子圖的控制權

for i = 1:length(targetChs1)
    % A. 拼接資料
    % 假設 data 為您的 261*21 timetable
    rawCells = data.(targetChs1{i})(startRow:endRow);
    combinedSig = vertcat(rawCells{:});
    
    % B. 精確裁切
    relativeStartIdx = round((startTime - (startRow-1)*secPerRow) * fs) + 1;
    segment = combinedSig(relativeStartIdx : relativeStartIdx + numSamples - 1);
    
    % C. 繪製子圖
    ax(i) = subplot(length(targetChs1), 1, i);
    plot(timeVec, segment, 'Color', [0.15 0.15 0.15]);
    
   % D. 格式調整
    ylabel(targetChs1{i}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    grid on;
    
    % 強制設定 Y 軸 (加上 'manual' 模式避免自動調整)
    set(gca, 'YLim', yLimit, 'YLimMode', 'manual'); 
    
    if i < length(targetChs1), set(gca, 'XTickLabel', []); end
end

xlabel('Time (seconds)');
sgtitle(['Raw EEG Data (No Filter) - Starting at ', num2str(startTime), 's']);


figure('Color', 'w', 'Name', 'Raw EEG Signal (954s-1014s)');
ax = zeros(1, length(targetChs2)); % 儲存每個子圖的控制權

for i = 1:length(targetChs2)
    % A. 拼接資料
    % 假設 data 為您的 261*21 timetable
    rawCells = data.(targetChs2{i})(startRow:endRow);
    combinedSig = vertcat(rawCells{:});
    
    % B. 精確裁切
    relativeStartIdx = round((startTime - (startRow-1)*secPerRow) * fs) + 1;
    segment = combinedSig(relativeStartIdx : relativeStartIdx + numSamples - 1);
    
    % C. 繪製子圖
    ax(i) = subplot(length(targetChs2), 1, i);
    plot(timeVec, segment, 'Color', [0.15 0.15 0.15]);
    
   % D. 格式調整
    ylabel(targetChs2{i}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    grid on;
    
    % 強制設定 Y 軸 (加上 'manual' 模式避免自動調整)
    set(gca, 'YLim', yLimit, 'YLimMode', 'manual'); 
    
    if i < length(targetChs2), set(gca, 'XTickLabel', []); end
end

xlabel('Time (seconds)');
sgtitle(['Raw EEG Data (No Filter) - Starting at ', num2str(startTime), 's']);
