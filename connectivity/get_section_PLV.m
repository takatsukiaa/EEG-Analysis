function plv_matrix = get_section_PLV(phase_matrix, plot)
%GET_SECTION_PLV Summary of this function goes here
%   Detailed explanation goes here
channel_names = fieldnames(phase_matrix);
num_ch = length(channel_names);
plv_matrix = ones(num_ch, num_ch); % 對角線必為 1
num_samples = length(phase_matrix.C3_AVG);

for i = 1:num_ch
    for j = i+1:num_ch
        p1 = phase_matrix.(channel_names{i});
        p2 = phase_matrix.(channel_names{j});
        
        e = exp(1i * (p1-p2));
        val = abs(sum(e)) / num_samples;
        plv_matrix(i, j) = val;
        plv_matrix(j, i) = val;
    end
end
mask = triu(true(size(plv_matrix)), 1);
plv_lower = plv_matrix;
plv_lower(mask) = NaN;
if plot
    % 使用 imagesc 畫出熱圖
    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [0, 0, 1200, 1200]); 
    
    % 2. 繪製熱圖
    h = heatmap(channel_names, channel_names, plv_lower);
    h.Title = 'Section PLV Values';
    h.Colormap = jet;
    h.ColorLimits = [0 1];
    h.CellLabelFormat = '%.2f';
    h.MissingDataColor = [1 1 1]; % Set NaN areas to white to hide them
    h.FontSize = 8; % 調整字體大小
    
    output_folder = 'PLV_Results_Images/No2_segment2_post_slice_5';
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    % 3. 高解析度存檔 (300 DPI)
    filename = fullfile(output_folder, 'PLV_Slice5_Segment2_post.png');
    exportgraphics(fig, filename, 'Resolution', 300);
    
    close(fig); % 關閉以釋放記憶體
end
    end

