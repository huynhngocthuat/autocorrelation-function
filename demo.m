clc;
clear;  
close all;
%Standard frequency
F0mean1 = 122;
F0std1 = 18;
F0mean2 = 217;
F0std2 = 23;
F0mean3 = 113;
F0std3 = 26;
F0mean4 = 232;
F0std4 = 40;

%Find frequency by autocorrelation function
[wav1, Fs1] = audioread('E:\Audacity\TinHieuHuanLuyen\phone_M1.wav');
[wav2, Fs2]  = audioread('E:\Audacity\TinHieuHuanLuyen\phone_F1.wav');
[wav3, Fs3] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_M1.wav');
[wav4, Fs4] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_F1.wav');
F1 = findFundamentalFrequency(wav1, Fs1);
F2 = findFundamentalFrequency(wav2, Fs2);
F3 = findFundamentalFrequency(wav3, Fs3);
F4 = findFundamentalFrequency(wav4, Fs4);
%Different
d1 = abs(F1 - F0mean1);
d2 = abs(F2 - F0mean2);
d3 = abs(F3 - F0mean3);
d4 = abs(F4 - F0mean4);
function f  = findFundamentalFrequency(x1, fs)
start = 30000;
x = x1(start:start + 1023); %Take 1024 samples
i = 1:1024;
t = i/fs;
y = xcorr(x); % Find autocorrelation 
y1 = y(1024:1024 + 600); %Take 600 samples autocorrelation

figure;
subplot(2,1,1);plot(t,x);
title('Audio Signal');
xlabel('Time');
ylabel('x');

M = movmean(y1,25); %Taking moving avarage of y1

subplot(2,1,2);plot(M);
title('Autocorrelation function');
xlabel('no of bins');
ylabel('M');
%Find maxima
n = 40;
global maxima;
global k;
k = 0;
while 1
    k = k + 1;
    while 1 
        if(M(n) < M(n-1))
            n = n +1;
        else
            m = n;
            break;
        end
    end

    while 1
        if(M(m) > M(m-1))
            m = m+1;
        else
            maxima(k) = m;
            break;
        end
    end
    if(m>400)
        break
    end
    n = m;
end
%Find max(maxima)
max = M(maxima(1));
index = maxima(1);
for h = 2 : 7
    if(M(maxima(h)) > max)
        max = M(maxima(h));
        index = maxima(h);
    end
end

f = fs/index;
end
