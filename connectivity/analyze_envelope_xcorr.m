function analyze_envelope_xcorr(eeg_struct, Fs, output_filename)
    % 取得所有通道名稱與數量
    channels = fieldnames(eeg_struct);
    num_ch = length(channels);
    
    % 若通道數不足2個，則無法配對
    if num_ch < 2
        error('結構體中的通道數量不足，無法進行配對分析！');
    end
    
    disp(['開始處理，共偵測到 ', num2str(num_ch), ' 個通道...']);
    
    % =========================================================
    % 階段 1：一次性計算所有通道的 Slow-gamma 能量包絡線 (極大提升效能)
    % =========================================================
    disp('階段 1: 正在萃取所有通道的 slow-gamma (25-50Hz) 包絡線...');
    
    window = round(Fs * 0.5);
    noverlap = round(Fs * 0.45);
    nfft = 4000;
    freq_range = [25, 50]; % 限定為 slow-gamma 頻段
    
    % 建立一個暫存的 struct 來放算好的包絡線
    env_data = struct(); 
    
    for i = 1:num_ch
        ch_name = channels{i};
        y = eeg_struct.(ch_name);
        
        % 時頻轉換與正規化
        [S, F, ~] = spectrogram(y, window, noverlap, nfft, Fs);
        S_norm = abs(S) ./ max(abs(S), [], 1);
        
        % 抓取 slow-gamma 頻段並平滑化
        idx = (F >= freq_range(1) & F <= freq_range(2));
        env = smoothdata(mean(S_norm(idx, :), 1), 'gaussian', 20);
        
        % 存入暫存結構體中
        env_data.(ch_name) = env;
    end
    
    % =========================================================
    % 階段 2：進行 C(N, 2) 兩兩配對互相關，並寫入 txt 檔
    % =========================================================
    disp(['階段 2: 正在計算兩兩互相關，並輸出至 ', output_filename, ' ...']);
    
    % 開啟 txt 檔案準備寫入 ('w' 代表覆寫模式)
    fid = fopen(output_filename, 'w');
    if fid == -1
        error('無法建立或開啟 txt 檔案，請檢查資料夾權限或路徑！');
    end
    
    % 寫入標題行 (使用 Tab 分隔，方便之後貼到 Excel)
    fprintf(fid, 'Channel_A\tChannel_B\tR_value\tLag_point\n');
    
    pair_count = 0;
    
    % ★ 核心邏輯：C(N, 2) 不重複配對迴圈
    for i = 1:(num_ch - 1)
        chA = channels{i};
        envA = env_data.(chA);
        
        % 內層迴圈從 i+1 開始，確保不重複且不和自己配對
        for j = (i + 1):num_ch
            chB = channels{j};
            envB = env_data.(chB);
            
            % 計算互相關與延遲
            [C, lag] = xcorr(envA, envB, 'coeff');
            [max_R, max_I] = max(C);
            max_lag = lag(max_I);
            
            % 寫入 txt 檔案 (保留 4 位小數)
            fprintf(fid, '%s\t%s\t%.4f\t%d\n', chA, chB, max_R, max_lag);
            
            pair_count = pair_count + 1;
        end
    end
    
    % 關閉檔案 (非常重要，否則檔案會被鎖住)
    fclose(fid);
    
    disp(['分析完成！共計產生 ', num2str(pair_count), ' 組配對結果。']);
end