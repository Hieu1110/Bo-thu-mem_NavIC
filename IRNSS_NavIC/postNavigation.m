function [navSolutions, eph,TOW] = postNavigation(trackResults, settings)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

if (settings.msToProcess < 50000)
    % Show the error message and exit
    disp('Record is to short or too few satellites tracked. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Find preamble start positions ==========================================

[subFrameStart, activeChnList] = findPreambles(trackResults, settings);

%% Decode ephemerides =====================================================

for channelNr = activeChnList

    %=== Convert tracking output to navigation bits =======================

    %--- Copy sub-frames long record from tracking output ---------------
    % -20 ??? 
    navBitsSamples = trackResults(channelNr).I_P(subFrameStart(channelNr) : ...
        subFrameStart(channelNr) + (600 * 4 * 20) - 1)';

    %--- Group every 20 vales of bits into columns ------------------------
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 1) / 20));

    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples);


    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    if navBits(1) > 0
        navBits = (navBits > 0);
    else
        navBits = (navBits < 0);
    end

    % disp(navBits);

    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.

    navBits = decodeFEC(navBits);

    % disp(navBits);

    crcCheck = checkCRC(navBits);
    if(crcCheck > 1)
        continue
    end

    navBitsBin = dec2bin(navBits);

    %=== Decode ephemerides and TOW of the first sub-frame ================
    [eph(trackResults(channelNr).PRN), TOW] = ...
        ephemeris(navBitsBin(1:292*4)');

    % disp(['Giá trị các tham số điều hướng của PRN ', int2str(trackResults(channelNr).PRN), ': ']);
    % disp(eph(trackResults(channelNr).PRN));

    %--- Exclude satellite if it does not have the necessary nav data -----
    if (isempty(eph(trackResults(channelNr).PRN).IODEC))
        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
    end
end

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Initialization =========================================================

% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.
readyChnList = activeChnList;

transmitTime = TOW;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################
prevClkDrift=0;
bias=0;
%Start Clock in samples
%Find the seconds nearest the startOfFrame and converting it to samples
%tmp=ceil(trackResults(1).absoluteSample(subFrameStart(1))/settings.samplingFreq);
%startClock_PVT_GPS=tmp*settings.samplingFreq;
%% Initialization of current measurement ==================================
for currMeasNr = 1%:fix((settings.msToProcess - max(subFrameStart)) / settings.navSolPeriod)

    % Exclude satellites, that are belove elevation mask
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
        readyChnList);
    
    if isempty(activeChnList)
        continue
    end
    % disp(['currMeasNr = ', int2str(currMeasNr)]); 

    % activeChnList = readyChnList;

    % Save list of satellites used for position calculation
    navSolutions.channel.PRN(activeChnList, currMeasNr) = ...
        [trackResults(activeChnList).PRN];

    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.channel.el(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    navSolutions.channel.az(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);


    %% THUAN
    if currMeasNr>1
        settings.startOffset=settings.startOffset-1000*(navSolutions.dt(currMeasNr-1))/settings.c;
    end

    navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
        trackResults, ...
        subFrameStart+settings.navSolPeriod*(currMeasNr-1), ...
        activeChnList, settings);

    % disp(['Giả khoảng cách từ bộ thu đến vệ tinh: ', num2str(navSolutions.channel.rawP(2, 1)), ' mét.'])

    xx=subFrameStart+settings.navSolPeriod*(currMeasNr-1);
    yy=inf(1,settings.numberOfChannels);
    if size(activeChnList,1)>1
        activeChnList=activeChnList';
    end
    for channelNr = activeChnList
        yy(channelNr) = (1023-trackResults(channelNr).remCodePhase(xx(channelNr)-1))/trackResults(channelNr).codeFreq(xx(channelNr)-1)*settings.c;
    end
    navSolutions.channel.rawP(:, currMeasNr)=navSolutions.channel.rawP(:, currMeasNr)+yy';

    %% END  THUAN

    %% Find satellites positions and clocks corrections =======================
    [satPositions, satClkCorr] = satpos(transmitTime+0, ...
        [trackResults(activeChnList).PRN], ...
        eph, settings);

    % disp("Vị trí của vệ tinh trên hệ tọa độ ECEF: ");
    % disp(['  X: ', num2str(satPositions(1,1))]);
    % disp(['  Y: ', num2str(satPositions(2,1))]);
    % disp(['  Z: ', num2str(satPositions(3,1))]);
    % disp("Thời gian tín hiệu truyền đến bộ thu: ");disp(satClkCorr(1));

    %% Find receiver position =================================================

    % 3D receiver position can be found only if signals from more than 3
    % satellites are available
    % if size(activeChnList, 2) > 3
    if numel(activeChnList) > 3

        %=== Calculate receiver position ==================================
        ps_fract1ms=rem(navSolutions.channel.rawP(activeChnList, currMeasNr)',settings.c*0.001);
        ps_1ms=floor(navSolutions.channel.rawP(activeChnList, currMeasNr)'/(settings.c*0.001));
        [xyzdt, ...
            navSolutions.channel.el(activeChnList, currMeasNr), ...
            navSolutions.channel.az(activeChnList, currMeasNr), ...
            navSolutions.DOP(:, currMeasNr)] = ...
            leastSquarePos(satPositions, ...
            ps_fract1ms,ps_1ms, satClkCorr, ...
            settings);

        %--- Save results -------------------------------------------------
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);
        navSolutions.dt(currMeasNr) = xyzdt(4);

        % Update the satellites elevations vector
        satElev = navSolutions.channel.el(:, currMeasNr);

        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.channel.correctedP(activeChnList, currMeasNr) = ...
            navSolutions.channel.rawP(activeChnList, currMeasNr) + ...
            satClkCorr' * settings.c + navSolutions.dt(currMeasNr);

        %% Coordinate conversion ==================================================

        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr), ...
            navSolutions.height(currMeasNr)] = cart2geo(...
            navSolutions.X(currMeasNr), ...
            navSolutions.Y(currMeasNr), ...
            navSolutions.Z(currMeasNr), ...
            5);


        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr));

        [navSolutions.E(currMeasNr), ...
            navSolutions.N(currMeasNr), ...
            navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
            xyzdt(3), ...
            navSolutions.utmZone);

    else % if size(activeChnList, 2) > 3
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
            ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;

        navSolutions.channel.az(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        navSolutions.channel.el(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));

        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.

    end % if size(activeChnList, 2) > 3

    %=== Update the transmit time ("measurement time") ====================
    transmitTime = transmitTime+ settings.navSolPeriod / 1000;

end %for currMeasNr...
