function frmError = checkCRC(navBits)


%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Detect errors in input data using CRC ---------------------------
[~,frmError] = step(crcDet,navBits');

