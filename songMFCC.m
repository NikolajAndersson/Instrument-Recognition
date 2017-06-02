%% Instrument classification

clear; clc; close all;


N = 2^13;
hop = N/2;
sr = 48000;

filepathsax = 'SMSAdata/sax/';
filepathvio = 'SMSAdata/violin/';
filepathcla = 'SMSAdata/clarinet/';
filepathtru = 'SMSAdata/trumpet/';

filename = ['0000'; '0001'; '0002'; '0003'; '0004'; '0005'; '0006'; '0007'; '0008'; '0009'];


%%
[filenames, pathname, filterindex] = uigetfile( '*.wav', 'WAV-files (*.wav)', 'Pick a file', 'MultiSelect', 'on');
%filenames = cellstr(filename);   %in case only one selected

%%
maxlim = 27000;
sampleM = [zeros(maxlim,length(filenames))];
for K = 1 : length(filenames)
  thisfullname = fullfile(pathname, filenames{K});
  [x,fs]=audioread(thisfullname);
  sampleM(:,K) = [x; zeros(maxlim-length(x),1)];
end
%%



% filter
T = triFilterBank(N, sr);
coef = 13; 
dataAmount = 10; 
hop = N/2;
instrumentAmount = 4; 
%%
saxdata = [];
saxID = 'sax';
figure; hold on; 
for i = 1:1
[s, ~] = audioread([filepathsax filename(i, :) '.wav']);
mfcc = getMFCCSong(s, N, T, coef);
saxdata(i,:) = [mfcc];
imagesc(saxdata(i,:))
end
title('Saxophone MFCC')