%% MFCC
clear; clc; close all; 
% high pass filter 
[s, fs] = audioread('o.wav');
N = 2^12;
han = hanning(N);
x = s(10001:10001+N-1,1);
xn1 = 0;
y = zeros(size(x));

for i = 1:N
    y(i) = x(i) - 0.97*xn1;
    xn1 = x(i); 
end
plot(x,'k'); hold on; plot(y,'r'); hold off; 
%% 
% hanning, fft and absolute value
X = abs(fft(y.*han, N));
% linear spacing between frequency bins
freqBins = linspace(0, fs/2, N/2);
figure
plot(freqBins,X(1:N/2))

%% 
coefAmount = 42;
f = 133; 
% mel-frequency spacing

melscale = zeros(1,coefAmount);
melscale(1) = f;
for i = 2:coefAmount
    if i <= 12        
        f = f + 66.66; % linear spacing
        melscale(i) = f;
    else        
        f = f * 1.0712; % log spacing
        melscale(i) = f;
    end    
end

%% 

plot(freqBins,abs(X(1:N/2))); hold on; plot(melscale,0,'k*')

%% 
% filterbank

heights = zeros(1, coefAmount);
heights(1) =  1/ (melscale(2) -  0);
for i = 2:coefAmount-1
    heights(i) = 1/ (melscale(i+1) -  melscale(i-1));
end
heights(coefAmount) = 1/ (0 -  melscale(i));

%% 

figure; 
plot(freqBins,abs(X(1:N/2))); hold on; plot(melscale,heights,'k')

hold off;
%% triangle filter bank
T = zeros(coefAmount,N/2);

for i = 1:coefAmount
    if i-1 < 1
        past = 0;
        current = melscale(i);
        next = melscale(i+1);
    elseif i+1 > coefAmount
        current = melscale(i);
        past = melscale(i-1);
        next = 0;
    else
        current = melscale(i);
        past = melscale(i-1);
        next = melscale(i+1);
    end
    
for j=1:N/2
    k = freqBins(j);
    
    if past <= k && k <= current % f(m-1) <= k, k <= f(m)
         % upgoing triangle            
        a = (k - past)/(current -  past);
        T(i,j) = a;
    elseif current <= k && k <= next % f(m) <= k <= f(m+1)
        % downgoing triangle
        a = (next - k)/(next -  current);
        T(i,j) = a;
    end
end
end
figure 
hold on
for i = 1:coefAmount
plot(T(i,:))

end
hold off;

%% Filter height might be wrong
Ts = zeros(size(T));

 for i = 1:coefAmount
    Ts(i,:) = T(i,:) /max(T(i,:) );
    Ts(i,:) = Ts(i,:)*heights(i) ;
 end
 
figure; hold on;
for i = 1:coefAmount
    plot(Ts(i,:))

end
hold off;
% seems better
T = Ts;

%% 
figure; hold on;
for i = 1:coefAmount
    plot(T(i,:))

end
hold off;
%% 
% AVERAGE 1/N sum T * abs X 

meanX = 1/N * (T*X(1:N/2)); 
%% log 
logX = log10(meanX);
%% discrete cosine transform
dctX = dct(logX)
figure
plot(dctX)

%% take the first 13

mfcc = dctX(1:13);
figure
plot(mfcc)

