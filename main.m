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

    %Read audio
[x1, Fs1] = audioread('E:\Audacity\TinHieuKiemThu\phone_M2.wav');
[x2, Fs2]  = audioread('E:\Audacity\TinHieuKiemThu\phone_F2.wav');
[x3, Fs3] = audioread('E:\Audacity\TinHieuKiemThu\studio_M2.wav');
[x4, Fs4] = audioread('E:\Audacity\TinHieuKiemThu\studio_F2.wav');

    %Vector F0
F1 = findFundamentalFrequency(x1, Fs1, 80);
F2 = findFundamentalFrequency(x2, Fs2, 156);
F3 = findFundamentalFrequency(x3, Fs3, 86);
F4 = findFundamentalFrequency(x4, Fs4, 31);

        %Flot F0
    %Plot phone_M2
figure('Name','Phone-M2','NumberTitle','off');
subplot(3,1,1);
plot((1:length(x1))/Fs1, x1);
title('Input signal phone-M2 ');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(3,1,2);
stem((1:length(F1))*0.025,F1,'.','linestyle','none');
ylim([0 450]);
title('F0 uses autocorrelaton function');
xlabel('Time(s)');
ylabel('Hz');
subplot(3,1,3);
stem((1:length(F1))*0.025,medianFilter(F1,7),'.','linestyle','none');
ylim([0 450]);
title('F0 uses medium filters');
xlabel('Time(s)');
ylabel('Hz');
    %Plot phone_F2
figure('Name','Phone-F2','NumberTitle','off');
subplot(3,1,1);
plot((1:length(x2))/Fs2, x2);
title('Input signal phone-F2 ');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(3,1,2);
stem((1:length(F2))*0.025,F2,'.','linestyle','none');
ylim([0 450]);
title('F0 uses autocorrelaton function');
xlabel('Time(s)');
ylabel('Hz');
subplot(3,1,3);
stem((1:length(F2))*0.025,medianFilter(F2,7),'.','linestyle','none');
ylim([0 450]);
title('F0 uses medium filters');
xlabel('Time(s)');
ylabel('Hz');
    %Plot studio_M2
figure('Name','Studio-M2','NumberTitle','off');
subplot(3,1,1);
plot((1:length(x3))/Fs3, x3);
title('Input signal studio-M2 ');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(3,1,2);
stem((1:length(F3))*0.025,F3,'.','linestyle','none');
ylim([0 450]);
title('F0 uses autocorrelaton function');
xlabel('Time(s)');
ylabel('Hz');
subplot(3,1,3);
stem((1:length(F3))*0.025,medianFilter(F3,7),'.','linestyle','none');
ylim([0 450]);
title('F0 uses medium filters');
xlabel('Time(s)');
ylabel('Hz');
    %Plot studio_F2
figure('Name','Studio-F2','NumberTitle','off');
subplot(3,1,1);
plot((1:length(x4))/Fs4, x4);
title('Input signal studio-F2 ');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(3,1,2);
stem((1:length(F4))*0.025,F4,'.','linestyle','none');
ylim([0 450]);
title('F0 uses autocorrelaton function');
xlabel('Time(s)');
ylabel('Hz');
subplot(3,1,3);
stem((1:length(F4))*0.025,medianFilter(F4,7),'.','linestyle','none');
ylim([0 450]);
title('F0 uses medium filters');
xlabel('Time(s)');
ylabel('Hz');

    %Fundamental frequency
% [F01, Fstd1] = F0Average(F1);
% [F02, Fstd2] = F0Average(F2);
% [F03, Fstd3] = F0Average(F3);
% [F04, Fstd4] = F0Average(F4);

    %Fundamental frequency use median filter
[F01, Fstd1] = F0Average(medianFilter(F1,5));
[F02, Fstd2] = F0Average(medianFilter(F2,5));
[F03, Fstd3] = F0Average(medianFilter(F3,5));
[F04, Fstd4] = F0Average(medianFilter(F4,5));

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
function f  = findFundamentalFrequency(x, fs, frameIndexPlot)
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
        if(maxPeak > 0.38)%xx[n]*0.38 threshold(voice, unvoice).
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