session = ["01" "06" "12"];
title_txt = ["pre","middle","post"];
% session = ["19" "20"];
% title_txt = ["1 month","3 month"];

fs = 250;
nfft =2^ceil(log(fs)/log(2));
window =floor(fs*0.5);
overlap = floor(window*0.8);

path = uigetdir('','請選擇受試者資料夾');
folder = dir(path);

cnt = 1;

for i = 1: length(session)
    path_f = [path '\' char(session(i)) '\'];
    data = dir([path_f '*_raw.txt']); %load all _raw.txt files in the folder
    Folds = {data.name}';
    
    baseline_data = dir([path '\baseline\' '*_raw.txt']);
    baseline_folds = {baseline_data.name}';
    
    % plot baseline
    fname =  strcat([path '\baseline\' char(baseline_folds(str2num(session(i))))] );
    try
       data = load(fname);
    catch
       error = error + 1;
    end
    [Mu_Band_power,alpha_power,CH1_f,CH2_f,artifact_index] = nft_as(data,fs,1,1);
    [S1,F1,T1,P1] = spectrogram(CH1_f,window,overlap,nfft,fs,'yaxis');
    [S2,F2,T2,P2] = spectrogram(CH2_f,window,overlap,nfft,fs,'yaxis');
    subplot(14,3,i);
    surf(T1,F1,P1,'edgecolor','none'); 
    axis xy; axis tight; 
    axis([1 length(CH1_f)/fs-1 0 30]);
    colormap jet; view(0,90);caxis([0 10]);
                
    subplot(14,3,i+3);
    surf(T2,F2,P2,'edgecolor','none'); axis tight;    
    axis([1 250*360-1 0 30]);
%     axis([1 length(CH2_f)/fs-1 0 30]);
    colormap jet;
    caxis([0 10]);
    
    %plot trial
    num = 1;
    for trail = 1:length(Folds)
        fname =  strcat([path_f char(Folds(trail))] );
        try
           data = load(fname);
        catch
           error = error + 1;
           continue;
        end
        
        [Mu_Band_power,alpha_power,CH1_f,CH2_f,artifact_index] = nft_as(data,fs,1,2);

        cutpoints = length(data(1,:))-360*250;
        CH1_f(1:cutpoints) = [];    
        CH2_f(1:cutpoints) = [];
        
        [S1,F1,T1,P1] = spectrogram(CH1_f,window,overlap,nfft,fs,'yaxis');
        [S2,F2,T2,P2] = spectrogram(CH2_f,window,overlap,nfft,fs,'yaxis');

        subplot(12,3,(trail-1)*6+1+(i-1));
        surf(T1,F1,P1,'edgecolor','none'); 
        axis xy; axis tight; 
        axis([1 length(CH1_f)/fs-1 0 30]);
        colormap jet; view(0,90);caxis([0 10]);
               
        if trail == 1        
            title(title_txt(i));
%             title(['session' char(session(i))]);
        end         
        
        subplot(12,3,(trail-1)*6+4+(i-1));
        surf(T2,F2,P2,'edgecolor','none'); 
        axis xy; axis tight; 
        axis([1 length(CH1_f)/fs-1 0 30]);
        colormap jet; view(0,90);caxis([0 10]);
        
        if trail == 3            
            ylabel('Frequency(Hz)');
        end
        
        if trail == 6            
            xlabel('Time (Seconds)'); 
        else           
            xticklabels([]);
        end
    end 
end
    set(gcf,'outerposition',get(0,'screensize')); 
    saveas(gcf, [path '.png'],'png');