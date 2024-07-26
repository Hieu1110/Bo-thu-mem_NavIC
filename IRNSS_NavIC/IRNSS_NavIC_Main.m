% -------------------------------------------------------------------------------------------
close all;fclose all;
addpath include;
global isNewRun; isNewRun = 1; % flag to set whether a function needs to rerun or not (old results in MAT files are used instead)
global sampletype; sampletype = 2; % 2: IQ, 1: I
global f_sampling; f_sampling = 2e6;    %81e6; % sampling frequency [Hz]
global nominalfreq; nominalfreq = 0*-13.550e6; % IF frequency [Hz]
global code_rate; code_rate =  1.023e6;
global code_length; code_length = 1023;

global samplesPDI; samplesPDI = ceil(f_sampling*4e-3); % 4 ms for GALILEO

global CNo_WINDOW; % [ms] window for C/No estimation
CNo_WINDOW =  1000; % [ms] 100 ms window
global FIFO_IP; FIFO_IP = zeros(1,CNo_WINDOW/4);
global FIFO_QP; FIFO_QP = zeros(1,CNo_WINDOW/4);
global Freq_sum; Freq_sum = 0;
global Old_Freq_sum; Old_Freq_sum = 0;

global settings; settings = NavIC_settings();

% str = 'D:\NavIC\Data\2705\navic20240527_9h+ (4).bin';
% fid=fopen(str,'rb');
% if (fid == -1)
%     disp('Không thể mở tệp dữ liệu cần thiết, thoát...')
%     return
% end
% 
% %% Acquisition
% fseek(fid, settings.skipNumberOfBytes, 'bof');
% tmp = fread(fid, 2*11*settings.samplesPerCode, settings.dataType)';
% data=tmp(1:2:end)+1i*tmp(2:2:end);
% acqResults = acquisitionIQ(data,settings);
% plotAcquisition(acqResults);
% disp(acqResults)
% 
% channel = test_preRun(acqResults,settings);
% showChannelStatus(channel,settings);
% 
% %% Tracking
% startTime = now;
% disp (['   Tracking bắt đầu lúc ', datestr(startTime)]);
% [trackResults, channel] = tracking_V0_IQ(fid, channel, settings);
% plotTracking(1:settings.numberOfChannels, trackResults, settings);


% mat_file = './2705/NavIC_TrackResults20240527_9h+ (4).mat';
% save (mat_file, "trackResults");


%% Calculate navigation solutions
TrackData = load('2705/NavIC_TrackResults20240527_9h.mat','-mat');
trackResults = TrackData.trackResults;

disp('   Tính toán kết quả định vị...');
navSolutions = postNavigation(trackResults, settings);

% In ra thông tin 
disp("Vị trí bộ thu: ");
disp(navSolutions);
disp("Trong đó, DOP là các giá trị:");
disp(['   GDOP: ', num2str(navSolutions.DOP(1,1))])
disp(['   PDOP: ', num2str(navSolutions.DOP(2,1))])
disp(['   HDOP: ', num2str(navSolutions.DOP(3,1))])
disp(['   VDOP: ', num2str(navSolutions.DOP(4,1))])
disp(['   TDOP: ', num2str(navSolutions.DOP(5,1))])


% disp('   Processing is complete for this data block');

%% Plot all results ===================================================
% disp ('   Ploting results...');

plotNavigation(navSolutions, settings);

mapPlot(navSolutions);

% disp('Post processing of the signal is over.');
