function [ T ] = triFilterBank(N,sr)
% Triangle filter bank
% Can be used for creating MFCC

lowestFrequency = 133.3333333;
linearFilters = 13;
linearSpacing = 66.66666666;
logSpacing = 1.0711703;
logFilters = 27;
% Keep this around for later....
totalFilters = linearFilters + logFilters;
% Now figure the band edges. Interesting frequencies are spaced
% by linearSpacing for a while, then go logarithmic. First figure
% all the interesting frequencies.
freqs = lowestFrequency + (0:linearFilters-1)*linearSpacing;
freqs(linearFilters+1:totalFilters+2) = ...
    freqs(linearFilters) * logSpacing.^(1:logFilters+2);
%%
% Lower, center, and upper band
% edges are all consequtive interesting frequencies.
lower = freqs(1:totalFilters);
center = freqs(2:totalFilters+1);
upper = freqs(3:totalFilters+2);
% We now want to combine FFT bins so that each filter has unit
% weight, assuming a triangular weighting function. First figure
% out the height of the triangle, then we can figure out each
% frequencies contribution
mfccFilterWeights = zeros(totalFilters,N);
triangleHeight = 2./(upper-lower);
fftFreqs = (0:N-1)/N*sr;
for chan=1:totalFilters
    mfccFilterWeights(chan,:) = ...
        (fftFreqs > lower(chan) & fftFreqs <= center(chan)).* ...
        triangleHeight(chan).*(fftFreqs-lower(chan))/(center(chan)-lower(chan)) + ...
        (fftFreqs > center(chan) & fftFreqs < upper(chan)).* ...
        triangleHeight(chan).*(upper(chan)-fftFreqs)/(upper(chan)-center(chan));
end

figure; hold on;
for i = 1:totalFilters
    plot(mfccFilterWeights(i,:));
    
end
title('Triangle Filter Bank');
hold off; 

T = mfccFilterWeights;
end

