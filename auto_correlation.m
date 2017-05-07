function [ ff ] = autocorrelation(segment,fs, lowerlimit, upperlimit )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tau = floor([1/upperlimit, 1/lowerlimit]* fs);
t = tau(1):tau(2); 
s = tau(2) - tau(1);
estimates = autocorr(segment,tau(2));
% 
% estimates = [];
% n = 1;
% for i = tau(1):tau(2)
%     delay = [zeros(i,1); segment];
%     estimates(:,n) = xcorr(segment(i+1:end), delay);
%     n = n + 1;
% end

[M, index] = max(estimates(tau(1):end));
%[M, index] = max(estimates);
% plot(t, estimates)
% hold on 
% plot(index,M,'o')

ff = 1/(t(index)/fs);