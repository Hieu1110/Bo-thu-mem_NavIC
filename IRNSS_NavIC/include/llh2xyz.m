% llh2xyz.m
%
% DESCRIPTION:
% Convert lat, long, height in WGS84 to ECEF X,Y,Z
% lat and long given in decimal degrees.
% altitude should be given in meters
%
% PFMS Project, 2009
% Nicholas Rutherford

% NOTES:
% None.

function [X,Y,Z] = llh2xyz(lat,long, h) 
  lat = lat/180*pi; %converting to radians
  long = long/180*pi; %converting to radians
  a = 6378137.0; % earth semimajor axis in meters 
  f = 1/298.257223563; % reciprocal flattening 
  e2 = 2*f -f^2; % eccentricity squared 
 
  chi = sqrt(1-e2*(sin(lat)).^2); 
  X = (a./chi +h).*cos(lat).*cos(long); 
  Y = (a./chi +h).*cos(lat).*sin(long); 
  Z = (a*(1-e2)./chi + h).*sin(lat);