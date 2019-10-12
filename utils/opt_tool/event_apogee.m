function [value,isterminal,direction] = event_apogee(t,Y,settings,varargin)
% Event function to stop simulation at apogee

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

Q = Y(10:13)';

%Inertial Frame velocities
vels = quatrotate(quatconj(Q),Y(4:6)');

%Stop checking if I'm in Propulsion Phase
if t > settings.tb
    value = vels(3);
else
    value = 1;
end

isterminal = 1;
direction = 1;


end

