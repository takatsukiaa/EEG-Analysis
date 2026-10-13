%%
clc;clear all;close all
%%
Fs=200;
windowlegth=200;
[FFT_y1,faxis]=pwelch(F3_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y2,faxis]=pwelch(F4_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y3,faxis]=pwelch(F7_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y4,faxis]=pwelch(F8_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y5,faxis]=pwelch(T3_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y6,faxis]=pwelch(T4_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y7,faxis]=pwelch(T5_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y8,faxis]=pwelch(T6_segment,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y9,faxis]=pwelch(i10,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y10,faxis]=pwelch(j10,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y11,faxis]=pwelch(k10,hamming(windowlegth),[],windowlegth,Fs);
[FFT_y12,faxis]=pwelch(l10,hamming(windowlegth),[],windowlegth,Fs);
%%
plot(faxis,FFT_y1)
%%
FFT=[FFT_y1(1:60) FFT_y2(1:60) FFT_y3(1:60) FFT_y4(1:60) FFT_y5(1:60) FFT_y6(1:60) FFT_y7(1:60) FFT_y8(1:60)];%FFT_y5(1:60)  FFT_y6(1:60) FFT_y3(1:60) FFT_y4(1:60) FFT_y5(1:60) FFT_y6(1:60) FFT_y7(1:51) FFT_y8(1:51)];