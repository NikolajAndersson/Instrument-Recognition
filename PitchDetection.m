%% Autocorrelation for seqment

[x, fs] = audioread('09viola.flac');

% Implement the autocorrelation pitch estimation method in, e.g., MATLAB as a function.
% The function should have the following input

l = 100;
h = 1000;
segment = x(fs:fs*2, 1);
autocorrelation(segment, fs, l, h)
%% Autocorrelation for entire signal
% frequencies in Hz
l = 10; % low limit
h = 2000; % high limit
segmentSize = 0.2; % in seconds
overlap = 100; % in percent
frequencies = pitchPlot('09viola.flac', segmentSize, overlap, l, h);
plot(frequencies)

%% 
algorithm = 'autocorrelation'

modelOrder = 10;
l = 10; % low limit
h = 2000; % high limit
segmentSize = 0.15; % in seconds
overlap = 10; % in percent
frequencies = pitchPlot2(algorithm,'09viola.flac', segmentSize, overlap, l, h, modelOrder);
plot(frequencies)
