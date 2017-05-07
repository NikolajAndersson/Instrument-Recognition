%% Autocorrelation for seqment
clc
clear 
[x, fs] = audioread('09viola.flac');

% Implement the autocorrelation pitch estimation method in, e.g., MATLAB as a function.
% The function should have the following input

l = 20;
h = 100;
winSize = 256
s = fs:fs+winSize-1;
segment = x(s, 1);
fundamentalfrequency(segment, fs, l, h)
%% Autocorrelation for entire signal
% frequencies in Hz
l = 10; % low limit
h = 2000; % high limit
segmentSize = 0.2; % in seconds
overlap = 100; % in percent
frequencies = pitchPlot('09viola.flac', segmentSize, overlap, l, h);
plot(frequencies)

%% 
l = 10; % low limit
h = 2000; % high limit
name = '09viola.flac';

segmentSize = 0.2; % in seconds
overlap = 100; % in percent
frequencies = pitchPlot(name, segmentSize, overlap, l, h);
plot(frequencies)
title(name, 'fontSize',16)
