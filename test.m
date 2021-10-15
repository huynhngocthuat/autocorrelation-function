clc;
clear;
close all;

a = [1,3,2,4,5,9,7,8,6];

% y = movmean(a,4);

[x, Fs1] = audioread('E:\Audacity\TinHieuHuanLuyen\phone_M1.wav');

% x = [1, 2, 1, 1];
y = xcorr(x);
t = 0.02;
samples = t * 16000; 
x1 = x(30000:30000 + 600);
y1 = autocorrelation(x1);
y2 = amdf(x1);
subplot(2,1,1);plot(y1);
title('Autocorrelation function');
xlabel('Index of sample');
ylabel('M');
subplot(2,1,2);plot(y2);
title('AMDF');
xlabel('Index of sample');
ylabel('M');





y2 = autocorrelation(x);
function y = amdf(x)
    y = zeros(1,length(x));
    frameLength = length(x);
    for lag = 0 : frameLength 
        for j = 1 : frameLength - lag
            y(lag+1) = y(lag+1) + abs(x(j) - x(j + lag));
        end                               
    end
end
function y = autocorrelation(x)
    y = zeros(1,length(x));
    frameIndex = 0;
    frameLength = length(x);
    for lag = 0 : frameLength 
        for j = 1 : frameLength - lag
            y(lag+1) = y(lag+1) + x(frameIndex + j) * x(frameIndex + j + lag);
        end                               
    end
end
function m = medianFilter(s, w)
    w2 = floor(w/2);
    w = 2*w2 + 1;

    n = length(s);
    m = zeros(w,n+w-1);
    s0 = s(1); sl = s(n);

    for i=0:(w-1)
        m(i+1,:) = [s0*ones(1,i) s sl*ones(1,w-i-1)];
    end
    m = median(m);
    m = m(w2+1:w2+n);
end
% figure;
% subplot(4,1,1);
% stem(F1, '.', 'LineStyle', 'none');
% xlabel('Frame');
% ylabel('Hz');
% subplot(4,1,2);
% stem(F2, '.', 'LineStyle', 'none');
% xlabel('Frame');
% ylabel('Hz');
% subplot(4,1,3);
% stem(F3, '.', 'LineStyle', 'none');
% xlabel('Frame');
% ylabel('Hz');
% subplot(4,1,4);
% stem(F4, '.', 'LineStyle', 'none');
% xlabel('Frame');
% ylabel('Hz');


