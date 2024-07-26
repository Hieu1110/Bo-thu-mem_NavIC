function [eph, TOW] = ephemeris(bits)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 2400 bits. The first element in
%the array must be the first bit of a subframe. The subframe ID of the
%first subframe in the array is not important.
%
%Function does not check parity!
%
%[eph, TOW] = ephemeris(bits)
%
%   Inputs:
%       bits        - bits of the navigation messages (5 subframes).
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%   Outputs:
%       TOW         - Time Of Week (TOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson

%% Check if the parameters are strings ====================================
if ~ischar(bits)
    error('The parameter BITS must be a character array!');
end


% Pi used in the GPS ~ NavIC coordinate system
gpsPi = 3.1415926535898; 
eph=struct('weekNumber',[],'deltan',[],'C_uc',[],'C_us',[],'C_ic',[],'C_is',[],'C_rc',[],'C_rs',[],'iDot',[],'M_0',[],...
    't_oe',[],'e',[],'sqrtA',[],'omega_0',[],'omega',[],'omegaDot',[],'i_0',[],'IODEC',[], 'T_GD',[],'a_f2',[],'a_f1',[],'a_f0',[],'t_oc',[]);
%% Decode all 4 sub-frames ================================================
count_subframe12 = 0; % Đếm số lần xuất hiện subframe 1 và 2. Đủ cả subframe 1 va 2 thì out
for i = 1:4

    if count_subframe12 == 2
        break;
    end
    %--- "Cut" one sub-frame's bits ---------------------------------------
    subframe = bits(292*(i-1)+1 : 292*i);

    %--- Decode the sub-frame id ------------------------------------------
    % For more details on sub-frame contents please refer to GPS IS.
    subframeID = bin2dec(subframe(28:29)) + 1; %

    %--- Decode sub-frame based on the sub-frames id ----------------------
    % The task is to select the necessary bits and convert them to decimal
    % numbers. For more details on sub-frame contents please refer to GPS
    % ICD (IS-GPS-200D).
    if subframeID == 1      %--- It is subframe 1 -------------------------------------   
            eph.weekNumber  = bin2dec(subframe(31:40)) + 1024;
            eph.deltan      = twosComp2dec(subframe(115:136)) * 2^(-41) * gpsPi;
            eph.C_uc        = twosComp2dec(subframe(157:171)) * 2^(-28);
            eph.C_us        = twosComp2dec(subframe(172:186)) * 2^(-28);
            eph.C_ic        = twosComp2dec(subframe(187:201)) * 2^(-28);
            eph.C_is        = twosComp2dec(subframe(202:216)) * 2^(-28);
            eph.C_rc        = twosComp2dec(subframe(217:231)) * 2^(-4);
            eph.C_rs        = twosComp2dec(subframe(232:246)) * 2^(-4);
            eph.iDot        = twosComp2dec(subframe(247:260)) * 2^(-43) * gpsPi;
            eph.IODEC       = bin2dec(subframe(137:144));
            eph.T_GD        = twosComp2dec(subframe(107:114)) * 2^(-31);
            eph.a_f0        = twosComp2dec(subframe(41:62)) * 2^(-31);
            eph.a_f1        = twosComp2dec(subframe(63:78)) * 2^(-43);
            eph.a_f2        = twosComp2dec(subframe(79:86)) * 2^(-55);
            eph.t_oc        = bin2dec(subframe(91:106)) * 2^4;

            count_subframe12 = count_subframe12 + 1;

    elseif subframeID == 2  %--- It is subframe 2 -------------------------------------
            eph.M_0         = twosComp2dec(subframe(31:62)) * 2^(-31) * gpsPi;
            eph.e           = bin2dec(subframe(79:110)) * 2^(-33);
            eph.sqrtA       = bin2dec(subframe(111:142)) * 2^(-19);
            eph.t_oe        = bin2dec(subframe(63:78)) * 2^4;
            eph.omega_0     = twosComp2dec(subframe(143:174))* 2^(-31) * gpsPi;
            eph.i_0         = twosComp2dec(subframe(229:260)) * 2^(-31) * gpsPi;
            eph.omega       = twosComp2dec(subframe(175:206))* 2^(-31) * gpsPi;
            eph.omegaDot    = twosComp2dec(subframe(207:228)) * 2^(-41) * gpsPi;

            count_subframe12 = count_subframe12 + 1;
            
    end 

end % for all 4 sub-frames ...

%% Compute the time of week (TOW) of the first sub-frames in the array ====
% Also correct the TOW. The transmitted TOW is actual TOW of the next
% subframe and we need the TOW of the first subframe in this data block
% (the variable subframe at this point contains bits of the last subframe). 
TOW = bin2dec(subframe(9:25)) * 12 - 12*4;
