function [ mfcc ] = getMFCC(s, N, T, coef )
file = [];
% filter of all the silent parts of the sound file
for i = 1:size(s,1)
    if abs(s(i,1)) > 0.001
        file = [file; s(i,1)];
    end
end

[l,m] = size(file);
hop = N/2;
frnop = floor(l/hop-1);
mfccMat = [];
% loop through sound and get MFCC
for i=1:frnop
    fr=file((i-1)*hop+(1:N));
    
    preEmphasized = filter([1 -.97], 1, fr);
    mX = abs(fft(preEmphasized.*hamming(N)));
    filteredX = T*mX;
    log10X = log10(filteredX);

    % DCT
    mfccMat(:,i) = dct(log10X);
end
% output average MFCC from the sound file
mfcc = mean(mfccMat,2);
mfcc = mfcc(2:coef+1);
plot(mfcc(1:13))

end

