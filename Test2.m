clc;
clear;  
close all;

F0mean1 = 122;
F0std1 = 18;
F0mean2 = 217;
F0std2 = 23;
F0mean3 = 113;
F0std3 = 26;
F0mean4 = 232;
F0std4 = 40;
%Find frequency by autocorrelation function
[x, fs] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_F1.wav');

frameTime = 0.03;
frameLength = frameTime * fs;
totalFrame = floor(length(x)/frameLength);
global start;
start = 1;
f = zeros(1,totalFrame);
for frame = 1 : totalFrame - 1
    x1 = x(start:start + frameLength); %Take frameLength samples
    y = zeros(1,length(x1)); %Samples autocorrelation 
    %Find autocorrelation
    frameIndex = 0;
    for lag = 0 : frameLength 
        for j = 1 : frameLength - lag
            y(lag+1) = y(lag+1) + x1(frameIndex + j) * x1(frameIndex + j + lag);
        end                               
    end
    
    %Plot 2 frames
    M = movmean(y,25); %Taking moving avarage of y1
    if(frame == 61)
        subplot(2,1,1);plot(M);
        title('Autocorrelation function');
        xlabel('Index of sample');
        ylabel('M');
    end
    if(frame == 62)
        subplot(2,1,2);plot(M);
        title('Autocorrelation function');
        xlabel('Index of sample');
        ylabel('M');
    end
    
    %Find maxima
    maxPeak = 0;
    index = 0;
    for i = floor(fs/450) : floor(fs/70)
        if (M(i) > M(i -1) && M(i) > M(i+1) && M(i) >= maxPeak)   
            maxPeak = M(i);                
            index = i;
        end
    end
    if(maxPeak > M(1)*0.3)
        f(frame) = fs/index;
    end
    
    %Go to the next frame
    start = start + frameLength;
end
F0_TTQ = medfilt1(f,7);
F01, F01std = F0Average(f, 30, 30);
F02, F02std = F0Average(F0_TTQ, 30, 30);

%Find average fundamemtal frequency
function [F0, F0std] = F0Average(f, first, last)
    count = 0 ; 
    F0 = 0 ;
    F0std = 0;
    for i = index:length(f) - 10
        if f(i) > 70 && f(i) < 450            
            F0 = F0 + f(i); 
            count = count + 1;           
        end
    end
    F0 = F0/count; %Average fundamental frequency
    
    count = 0;
    for j = index:length(f) - 10
        if f(j) > 70 && f(j) < 450            
            F0std = F0std + abs(F0-f(j));
            count = count + 1;
        end
    end
    F0std = F0std/count; %Standard deviation
end 
