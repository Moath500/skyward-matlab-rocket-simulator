function [ uw,vw,ww ] = windgen( AzMin,AzMax,ElMin,ElMax,MagMin,MagMax )
%windgen( AzMin,AzMax,ElMin,ElMax,MagMin,MagMax )
%function that generates wind components in NED axes based on altitude
%
%Very basic implementation: it has to be improved (Power Spectral Density
%etc..)
%
%
%Vector Orientation
%AzMin = 0; Minimum angle of Azimuth from North
%AzMax = 2*pi; Maximum angle of Azimuth from North
%ElMin = 0; Minimum angle of Elevation
%ElMax = pi/2; Maximum angle of Elevatiom

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

%Generating random values for orientation and magnitude
Az = AzMin + (AzMax-AzMin)*rand;
El = ElMin + (ElMax-ElMin)*rand;
Mag = MagMin + (MagMax-MagMin)*rand;

%Random Wind Vector

%uw = Mag*cos(Az)*cos(El);
%vw = Mag*cos(El)*sin(Az);
%ww = -Mag*sin(El);

R = Mag*angle2dcm(Az,El,0,'ZYX');

uw = R(1,1);
vw = R(1,2);
ww = R(1,3);






end

