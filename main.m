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

    %Read audio
[x1, Fs1] = audioread('E:\Audacity\TinHieuHuanLuyen\phone_M1.wav');
[x2, Fs2]  = audioread('E:\Audacity\TinHieuHuanLuyen\phone_F1.wav');
[x3, Fs3] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_M1.wav');
[x4, Fs4] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_F1.wav');

    %Vector F0
F1 = findFundamentalFrequency(x1, Fs1, 115, 0.03);
F2 = findFundamentalFrequency(x2, Fs2, 18, 0.03);
F3 = findFundamentalFrequency(x3, Fs3, 70, 0.03);
F4 = findFundamentalFrequency(x4, Fs4, 96, 0.02);

    %Flot F0
figure;
subplot(2,1,1);
stem(F1, '.', 'LineStyle', 'none');
title('phone-M1');
xlabel('Frame');
ylabel('Hz');
subplot(2,1,2);
stem(F2, '.', 'LineStyle', 'none');
title('phone-F1');
xlabel('Frame');
ylabel('Hz');
figure;
subplot(2,1,1);
stem(F3, '.', 'LineStyle', 'none');
title('studio-M1');
xlabel('Frame');
ylabel('Hz');
subplot(2,1,2);
stem(F4, '.', 'LineStyle', 'none');
title('studio-F1');
xlabel('Frame');
ylabel('Hz');

    %Fundamental frequency
% [F01, Fstd1] = F0Average(F1, 20);
% [F02, Fstd2] = F0Average(F2, 10);
% [F03, Fstd3] = F0Average(F3, 20);
% [F04, Fstd4] = F0Average(F4, 10);

    %Fundamental frequency use median filter
[F01, Fstd1] = F0Average(medianFilter(F1,7), 20);
[F02, Fstd2] = F0Average(medianFilter(F2,7), 10);
[F03, Fstd3] = F0Average(medianFilter(F3,7), 20);
[F04, Fstd4] = F0Average(medianFilter(F4,7), 10);

    %Different
d1 = abs(F01 - F0mean1);
d2 = abs(F02 - F0mean2);
d3 = abs(F03 - F0mean3);
d4 = abs(F04 - F0mean4);
dstd1 = abs(Fstd1 - F0std1);
dstd2 = abs(Fstd2 - F0std2);
dstd3 = abs(Fstd3 - F0std3);
dstd4 = abs(Fstd4 - F0std4);
    %Find fundamental frequency by autocorrelation function
    %frameIndexPhot: plot frame khong tuan hoan
    %frameIndexPhot+1: plot frame tuan hoan
    %frameTime: xac dinh do dai cua frame
function f  = findFundamentalFrequency(x, fs, frameIndexPlot, frameTime)
    frameLength = floor(frameTime * fs);
    totalFrame = floor(length(x)/frameLength);
    global start;
    start = 1;
    f = zeros(1,totalFrame);
    
    for frame = 1 : totalFrame - 1
        x1 = x(start:start + frameLength); %Take frameLength samples
        y = zeros(1,length(x1)); %Samples autocorrelation 
        
        %Find autocorrelation
        for lag = 0 : frameLength 
            for j = 1 : frameLength - lag
                y(lag+1) = y(lag+1) + x1(j) * x1(j + lag);
            end                               
        end
    
        %Plot 2 frames
        M = movmean(y,25); 
        %Normalize
        for i = 2 : length(M)
            M(i) = M(i) / M(1);
        end
        M(1) = 1;
        if(frame == frameIndexPlot)
            figure;
            subplot(2,1,1);plot(M);
            title('non-cyclic');
            xlabel('Index of sample');
            ylabel('M');
        end
        if(frame == frameIndexPlot + 1)
            subplot(2,1,2);plot(M);
            title('cyclic');
            xlabel('Index of sample');
            ylabel('M');
        end

        %Find maxima
        maxPeak = 0;% globally maximal
        index = 0;
        for i = floor(fs/450) : floor(fs/70) % 70Hz -> 450Hz 
            if (M(i) > M(i -1) && M(i) > M(i+1) && M(i) >= maxPeak)   
                maxPeak = M(i);                
                index = i;
            end
        end
        if(maxPeak > M(1)*0.3)%xx[n]*0.3 threshold(voice, unvoice).
            f(frame) = fs/index;
        end

        %Go to the next frame
        start = start + frameLength;
    end

end

    %Function find standard deviation(F0std) and average fundamental frequency(F0)
function [F0, F0std] = F0Average(f, index)
    count = 0 ; 
    F0 = 0 ;
    F0std = 0;
    for i = index:length(f)-20
        if f(i) > 70 && f(i) < 450            
            F0 = F0 + f(i); 
            count = count + 1;           
        end
    end
    F0 = F0/count; %Average fundamental frequency
    
    count = 0;
    for j = index:length(f)-20
        if f(j) > 70 && f(j) < 450            
            F0std = F0std + abs(F0-f(j));
            count = count + 1;
        end
    end
    F0std = F0std/count; %Standard deviation
end 
    % Median filter
    %s: F0[]
    %w: Filter size
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