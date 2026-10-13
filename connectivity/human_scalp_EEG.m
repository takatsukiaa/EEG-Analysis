clc;
clear all;
close all;

filepath = 'No.1 Sz 15min54.03sec.edf'
%% 打開 EDF 檔
edfFileAnalyzer

%% 看 Annotation (需要edf+的檔案格式)
[data, anns] = edfread(filepath)
%% 算Sample rate

info = edfinfo(filepath);
Fs = info.NumSamples ./ seconds(info.DataRecordDuration);

