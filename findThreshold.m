clc;
clear;  
close all;

%Standard frequency
F0mean1 = 129;
F0std1 = 18.6;
F0mean2 = 145;
F0std2 = 33.7;
F0mean3 = 155;
F0std3 = 30.8;
F0mean4 = 200;
F0std4 = 46.1;
deviation = zeros(50,4);
standardDeviation = zeros(50,4);
count = 0;
for threshold = 0.01:0.01:0.5
    count = count + 1;
        %Read audio
    [x1, Fs1] = audioread('E:\Audacity\TinHieuHuanLuyen\phone_M1.wav');
    [x2, Fs2]  = audioread('E:\Audacity\TinHieuHuanLuyen\phone_F1.wav');
    [x3, Fs3] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_M1.wav');
    [x4, Fs4] = audioread('E:\Audacity\TinHieuHuanLuyen\studio_F1.wav');

        %Vector F0
    F1 = findFundamentalFrequency(x1, Fs1, threshold);
    F2 = findFundamentalFrequency(x2, Fs2, threshold);
    F3 = findFundamentalFrequency(x3, Fs3, threshold);
    F4 = findFundamentalFrequency(x4, Fs4, threshold);
        %Fundamental frequency
    [F01, Fstd1] = F0Average(F1);
    [F02, Fstd2] = F0Average(F2);
    [F03, Fstd3] = F0Average(F3);
    [F04, Fstd4] = F0Average(F4);

        %Different
    deviation(count,1) = abs(F01 - F0mean1);
    deviation(count,2) = abs(F02 - F0mean2);
    deviation(count,3) = abs(F03 - F0mean3);
    deviation(count,4) = abs(F04 - F0mean4);
    
    standardDeviation(count,1) = abs(Fstd1 - F0std1);
    standardDeviation(count,2) = abs(Fstd2 - F0std2);
    standardDeviation(count,3) = abs(Fstd3 - F0std3);
    standardDeviation(count,4) = abs(Fstd4 - F0std4);
end
minDeviation = deviation(1,1) + deviation(1,2) + deviation(1,3) + deviation(1,4);
minStandardDeviation = standardDeviation(1,1) + standardDeviation(1,2) + standardDeviation(1,3) + standardDeviation(1,4);
index = 1;
indexS = 1;
for i = 2 : 50
    sum = deviation(i,1) + deviation(i,2) + deviation(i,3) + deviation(i,4);
    sumStd = standardDeviation(i,1) + standardDeviation(i,2) + standardDeviation(i,3) + standardDeviation(i,4);
    if(sum < minDeviation)
        minDeviation = sum;
        index = i;
    end
    if(sumStd < minStandardDeviation)
        minStandardDeviation = sumStd;
        indexS = i;
    end
end
threshold = index / 100;
thresholdStd = indexS / 100;
figure;
subplot(2,1,1);plot(deviation);
title('Deviation');
xlabel('Threshold(%)');
ylabel('Delta');
subplot(2,1,2);plot(standardDeviation);
title('Standard Deviation');
xlabel('Threshold(%)');
ylabel('Delta');

    %Find fundamental frequency by autocorrelation function
    %frameIndexPhot: plot frame khong tuan hoan
    %frameIndexPhot+1: plot frame tuan hoan
function f  = findFundamentalFrequency(x, fs, threshold)
    frameTime = 0.025;
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

        %Find maxima
        maxPeak = 0;% globally maximal
        index = 0;
        for i = floor(fs/450) : floor(fs/70) % 70Hz -> 450Hz 
            if (M(i) > M(i -1) && M(i) > M(i+1) && M(i) >= maxPeak)   
                maxPeak = M(i);                
                index = i;
            end
        end
        if(maxPeak > M(1)*threshold)%xx[n]*0.3 threshold(voice, unvoice).
            f(frame) = fs/index;
        end

        %Go to the next frame
        start = start + frameLength;
    end

end

    %Function find standard deviation(F0std) and average fundamental frequency(F0)
function [F0, F0std] = F0Average(f)
    count = 0 ; 
    F0 = 0 ;
    F0std = 0;
    for i = 1:length(f)
        if f(i) > 70 && f(i) < 450            
            F0 = F0 + f(i); 
            count = count + 1;           
        end
    end
    F0 = F0/count; %Average fundamental frequency
    
    count = 0;
    for j = 1:length(f)
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