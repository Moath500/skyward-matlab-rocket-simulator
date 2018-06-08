function [uw,vw,ww] = wind_matlab_generator(settings,z,t,Hour,Day)
% wind_generator(settings,z,t)
% Function that generates wind components in NED reference frame
% Based on hwm07 model


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

if nargin == 3
    if settings.wind.HourMin == settings.wind.HourMax && settings.wind.HourMin == settings.wind.HourMax
    Day = settings.wind.DayMin;
    Hour = settings.wind.HourMin;
    end
end

Seconds = Hour*3600;

%% HORIZONTAL WIND

[uw,vw] = atmoshwm(settings.wind.Lat,settings.wind.Long,h,'day',Day,...
    'seconds',Seconds+t,'model','quiet','version','14');    %NED reference
ww = settings.wind.ww;


end


