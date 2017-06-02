function [ mfcc ] = getMFCC(s, N, T, coef )

% [l,m] = size(s);
% hopsize = N/2;
% 
% frnop = floor(l/hopsize-1);
% file = s(:,1);
% for i=1:frnop
%     fr=s((i-1)*hop+(1:N));
%     
%     preEmphasized = filter([1 -.97], 1, fr);
%     mX(:,i)=abs(fft(fr.*hanning(N))))';
% end

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
mfcc = mfcc(2:coef+1);
plot(mfcc(1:13))



end

