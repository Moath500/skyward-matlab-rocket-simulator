function [uw,vw,ww] = wind_generator(settings,z,t,Q)
% wind_generator( settings,z,t,Q )
% Function that generates wind components in body reference frame
% Based on hwm07 model
% nargin = 4:rotation of wind components in body axis
% nargin = 3: wind components in NED


% Author: Gabriele Poiana
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: gabriele.poiana@skywarder.eu
% Website: http://www.skywarder.eu
% January 2016; Last revision: 17.I.2016
% License:  2-clause BSD


h = -z+settings.z0;
if h < 0
    h = 0;
end

%% HORIZONTAL WIND

[uw,vw] = atmoshwm(settings.wind.Lat,settings.wind.Long,h,'day',settings.wind.Day,...
    'seconds',settings.wind.Seconds+t,'model','quiet','version','14'); %NED reference
ww = settings.wind.ww;

% Rotation in body reference frame (ascend)
if nargin == 4

    wind=quatrotate(Q,[uw vw ww]); % wind speed in body reference

    uw = wind(1);
    vw = wind(2);
    ww = wind(3);
end

end


