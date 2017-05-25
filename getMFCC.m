function [ mfcc ] = getMFCC(s, N, T, coef )

%N = 2^12;
%[s, sr] = audioread(name);
sample = s(40001:40001+N-1,1);

preEmphasized = filter([1 -.97], 1, sample);

% with hanning
mX2=abs(fft(preEmphasized.*hanning(N)));

 % filter + log10
filteredX =T*mX2;
log10X = log10(filteredX);
%plot(log10X)

% DCT
mfcc = dct(log10X);
mfcc = mfcc(1:coef);
%plot(mfcc(1:13))



end

