function [ T ] = triFilterBank(N,sr)
% Triangle filter bank with Mel-frequency spacing

lowestFrequency = 133.3333333;
linearFilters = 13;
linearSpacing = 66.66666666;
logSpacing = 1.0711703;
logFilters = 27;

% Total amount of filters 
totalFilters = linearFilters + logFilters;
% Find the filter frequencies
freqs = lowestFrequency + (0:linearFilters-1)*linearSpacing;
freqs(linearFilters+1:totalFilters+2) = ...
    freqs(linearFilters) * logSpacing.^(1:logFilters+2);

% Lower, center, and upper band for triangles
lower = freqs(1:totalFilters);
center = freqs(2:totalFilters+1);
upper = freqs(3:totalFilters+2);
% Step 5
mfccFilterWeights = zeros(totalFilters,N);
triangleHeight = 2./(upper-lower);
% N bins to frequencies
fftFreqs = (0:N-1)/N*sr;
% Calculate filters using the fft frequencies and the mel-frequencies
for i = 1:totalFilters
    mfccFilterWeights(i,:) = ...
        (fftFreqs > lower(i) & fftFreqs <= center(i)).* ...
        triangleHeight(i).*(fftFreqs-lower(i))/(center(i)-lower(i)) + ...
        (fftFreqs > center(i) & fftFreqs < upper(i)).* ...
        triangleHeight(i).*(upper(i)-fftFreqs)/(upper(i)-center(i));
end

figure; hold on;
%pRes = 1400;
for i = 1:totalFilters
    plot(fftFreqs, mfccFilterWeights(i,:));    
end
title('Triangle Filter Bank', 'FontSize',18);
ylabel('Magnitude','FontSize',14)
xlabel('Hz','FontSize',14)
hold off; 

T = mfccFilterWeights;
end

