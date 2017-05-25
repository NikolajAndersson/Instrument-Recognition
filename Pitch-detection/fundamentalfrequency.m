function [ ff ] = fundamentalfrequency(segment,fs, lowerlimit, upperlimit )
%Pitch detection of periodic signals. Only works for melodies, find the fundamental frequency. 
%   The algorithm uses auto correlation and outputs a frequency in Hz.

tau = floor([1/upperlimit, 1/lowerlimit]* fs);
t = tau(1):tau(2); 
s = tau(2) - tau(1);
estimates = autocorr(segment,tau(2));

[M, index] = max(estimates(tau(1):end));

ff = 1/(t(index)/fs);