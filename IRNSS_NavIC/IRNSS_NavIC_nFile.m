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
k = 1;
for i =  1:1:24     
    close all;fclose all;

    disp(i);

    % str = "D:\NavIC\Data\2705\navic20240527_9h+ (" + i + ").bin";
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
    % % disp(acqResults)
    % 
    % channel = test_preRun(acqResults,settings);
    % showChannelStatus(channel,settings);
    % 
    % %% Tracking
    % startTime = now;
    % disp (['   Tracking bắt đầu lúc ', datestr(startTime)]);
    % [trackResults, channel] = tracking_V0_IQ(fid, channel, settings);
    % plotTracking(1:settings.numberOfChannels, trackResults, settings);

    % mat_file = "./2705/NavIC_TrackResults20240527_9h+ (" + i + ")";
    % save (mat_file, "trackResults");


    %% Calculate navigation solutions
    mat_file = "2705/NavIC_TrackResults20240527_9h+ (" + i + ").mat";
    TrackData = load(mat_file);
    trackResults = TrackData.trackResults;

    disp('   Tính toán kết quả định vị...');
    navSolutions = postNavigation(trackResults, settings);

    % plotNavigation(navSolutions, settings);

    if isempty(navSolutions)
        continue
    else

        % Ensure the number of elements is 14 by padding with NaNs if necessary
        num_elements = 14;

        % Padding PRN array
        PRN_size = size(navSolutions.channel.PRN(:,1), 1);
        PRN_padding = num_elements - PRN_size;
        PRN = [navSolutions.channel.PRN(:,1); zeros(PRN_padding, 1)];

        % Padding el array
        el_size = size(navSolutions.channel.el(:,1), 1);
        el_padding = num_elements - el_size;
        el = [navSolutions.channel.el(:,1); NaN(el_padding, 1)];

        % Padding az array
        az_size = size(navSolutions.channel.az(:,1), 1);
        az_padding = num_elements - az_size;
        az = [navSolutions.channel.az(:,1); NaN(az_padding, 1)];

        % Initialize arrays for the struct with NaN
        sortedPRN = 1:14;
        sortedEl = NaN(1, num_elements);
        sortedAz = NaN(1, num_elements);

        for j = 1:numel(PRN)
            if PRN(j) > 0 && PRN(j) < num_elements
                sortedPRN(PRN(j)) = PRN(j);
                sortedEl(PRN(j)) = el(j);
                sortedAz(PRN(j)) = az(j);
            end
        end

        allNavSolutions.channel.PRN(:, k) = sortedPRN(:);
        allNavSolutions.X(k)           = navSolutions.X;
        allNavSolutions.Y(k)           = navSolutions.Y;
        allNavSolutions.Z(k)           = navSolutions.Z;
        allNavSolutions.dt(k)          = navSolutions.dt;
        allNavSolutions.DOP(:, k)      = navSolutions.DOP(:, 1);
        allNavSolutions.latitude(k)    = navSolutions.latitude;
        allNavSolutions.longitude(k)   = navSolutions.longitude;
        allNavSolutions.height(k)      = navSolutions.height;
        allNavSolutions.E(k)           = navSolutions.E;
        allNavSolutions.N(k)           = navSolutions.N;
        allNavSolutions.U(k)           = navSolutions.U;

        allNavSolutions.channel.az(:, k) = sortedAz(:);
        allNavSolutions.channel.el(:, k) = sortedEl(:);

        k = k+1;
    end


end

if ~isempty(allNavSolutions)
    plotNavigation(allNavSolutions, settings);
else
    disp("Không đủ dữ liệu!")
    return
end

% SkyPlot
figure;
hAxes = axes;

skyPlot1(hAxes, ...
    allNavSolutions.channel.az, ...
    allNavSolutions.channel.el, ...
    allNavSolutions.channel.PRN(:, 1));

title(hAxes, ['Sky plot (mean PDOP: ', ...
    num2str(mean(allNavSolutions.DOP(2,:))), ')']);

mapPlot(allNavSolutions);

