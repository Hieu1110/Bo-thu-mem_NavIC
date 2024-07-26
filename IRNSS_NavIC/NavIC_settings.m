function settings = NavIC_settings()
settings.msToProcess        = 60*1e3;        %[ms]   %%

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipNumberOfBytes  = 10*4*2e6;

settings.Brej = 5e3;
settings.relativeFreq   = 0;
settings.numberOfChannels   = 14;

settings.IF             = 0*-13.550e6;     %-13.550e6;      %[Hz]
settings.samplingFreq   = 2e6; %8184000;%    %[Hz]
settings.codeFreqBasis  = 1.023e6;      %[Hz]

settings.codeLength     = 1023;
settings.samplesPDI     = settings.samplingFreq*10e-3;

settings.samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
settings.dataType       = 'int16';
if strcmp(settings.dataType,'int8')
    settings.dataTypeSize   =  1;
else
    settings.dataTypeSize   =  2;
end

settings.skipAcquisition    = 2e6;

settings.acqSatelliteList   = 1:14;         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 16;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 1.4; % 2.4;

%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       %[Hz] 2
settings.dllCorrelatorSpacing    = 0.5;     %[chips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 25;      %[Hz] 25

%% Navigation solution settings ===========================================
% Period for calculating pseudoranges and position
settings.navSolPeriod   = 100;      %[ms]
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask  = 10;       %[degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr    = 1;        % 0 - Off
                                    % 1 - On
% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E = nan;
settings.truePosition.N = nan;
settings.truePosition.U = nan;


%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On
settings.plotAcquisition    = 1;

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 0; %68.802;       %[ms] Initial sign. travel time  ???

%% CNo Settings============================================================
% Accumulation interval in Tracking (in Sec)
settings.CNo.accTime = 0.004;
% Accumulation interval for computing VSM C/No (in ms)
settings.CNo.VSMinterval = 400;