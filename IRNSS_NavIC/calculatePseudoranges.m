function [pseudoranges] = calculatePseudoranges(trackResults, ...
                                                msOfTheSignal, ...
                                                channelList, settings)
%calculatePseudoranges finds relative pseudoranges for all satellites
%listed in CHANNELLIST at the specified millisecond of the processed
%signal. The pseudoranges contain unknown receiver clock offset. It can be
%found by the least squares position search procedure. 
%
%[pseudoranges] = calculatePseudoranges(trackResults, msOfTheSignal, ...
%                                       channelList, settings)
%
%   Inputs:
%       trackResults    - output from the tracking function
%       msOfTheSignal   - pseudorange measurement point (millisecond) in
%                       the trackResults structure
%       channelList     - list of channels to be processed
%       settings        - receiver settings
%
%   Outputs:
%       pseudoranges    - relative pseudoranges to the satellites. 

%--- Set initial travel time to infinity ----------------------------------
% Later in the code a shortest pseudorange will be selected. Therefore
% pseudoranges from non-tracking channels must be the longest - e.g.
% infinite. 
travelTime = inf(1, settings.numberOfChannels);

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

%--- For all channels in the list ... 
if(size(channelList,1)>1)
    channelList=channelList';
end


for channelNr = channelList

    sampleIndex = msOfTheSignal(channelNr); % Lấy chỉ số mẫu từ msOfTheSignal
    absoluteSampleValue = trackResults(channelNr).absoluteSample(sampleIndex); % Lấy giá trị mẫu tương ứng
    
    %--- Compute the travel times -----------------------------------------    
    travelTime(channelNr) = absoluteSampleValue / samplesPerCode;
end 
   

%--- Truncate the travelTime and compute pseudoranges ---------------------
minimum         = floor(min(travelTime));
travelTime      = travelTime - minimum + settings.startOffset;

%--- Convert travel time to a distance ------------------------------------
% The speed of light must be converted from meters per second to meters
% per millisecond. 
pseudoranges    = travelTime * (settings.c / 1000);

