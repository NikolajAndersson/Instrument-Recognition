%% Onset Detection
% import
[s, fs] = audioread('038_phrase_disco_simple_slow_sticks_ds.wav');
fid = fopen('onset.txt');
data = textscan(fid,'%f%s','delimiter',' ');
fclose(fid);

allOnsets = data{1};
instruments = data{2};
%% calc

step = 1/fs; 
t = 0:step:length(s)/fs - step;

bassdrum = [];
hihat = [];
snare = [];

for i = 1:length(allOnsets)
    if(strcmp(instruments(i), 'Bass'))
        bassdrum = [bassdrum, allOnsets(i)];
    end
    if(strcmp(instruments(i), 'Hihat'))
        hihat = [hihat, allOnsets(i)];
    end
    if(strcmp(instruments(i), 'Snare'))
        snare = [snare, allOnsets(i)];
    end
    
end


%% Plot
plot(t,s)

hold on
%plot(allOnsets,0, '*')
plot(bassdrum,0, 'r*')
plot(hihat,0, 'g*')
plot(snare,0, 'k*')

%% 
% choose a window function, hanning

% window length in samples, power of 2, corresponds to 35 ms
ms = 0.009;
win = 2^nextpow2(fs * ms); 
% choose hopsize as half of the window size. 
hop = win/2;

max_frames = floor(length(s)/hop - 1);
hopsInSec = (hop*(1:max_frames-1))/fs;

k = [0:win/2-1,  win/2:-1:1];
for i = 1:max_frames-1
    fr = s((i-1)*hop+(1:win)) .* hanning(win);
    X(:,i) = fft(fr,win);
    
    E(i) =  k * abs(X(:,i)).^2  / win;
    
end

%imagesc(abs(X))
%figure
En=E/max(abs(E-mean(E)));
%plot(t,s, 'k')
hold on;
plot(hopsInSec,En,'r')

%% 
for i=1:fn,
    e(i)=sum(R(:,i).^2)/w; % (4)
end
[eacf, elags]=xcorr(e);
figure
plot(elags*h/sr,eacf)