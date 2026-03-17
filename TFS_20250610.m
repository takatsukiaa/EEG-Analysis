%%
clc;clear all;close all% clear data and code
%%
%Time Frequency Spectrogram
Fs=200;% make sure your sampling rate
y= F3_segment12;% input your data name, and make sure the amplitude unit
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);% notch filter
y1 = filtfilt(d,y);% data filter 60Hz
%y2 = bandpass(y1,[25 50],Fs);
%y2 = bandstop(y1,[119 121],Fs);
[S,F,T]=spectrogram(y1,Fs*0.5,Fs*0.45,Fs*2,Fs);% Short-time Fourier Transform
%-----Normalize, normalize to the max value of the time
[aa,bb]=size(S);
Smax=max(S);
for k=1:bb
    for kk=1:aa
        S(kk,k)=S(kk,k)/Smax(k);
    end
end
%%
% ----plot figure
figure
surf(T,F,abs(S),'EdgeColor','none');
colormap jet
colorbar
set(gca,'Clim',[0,1],'FontSize',15,'FontWeight','bold')% color bar axial set
az=0;
el=90;
view(az,el);
ylim([20 50]);
xlim([0 inf]);
ylabel('Frequency (Hz)')
xlabel('time(s)')
set(gca,'FontSize',16,'FontWeight','bold','linewidth',1);
box off
