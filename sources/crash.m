function [ value,isterminal,direction ] = crash( t,Y,settings,varargin )
%Event Function for ODE to stop simulation @ landing

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD
value = Y(3);
isterminal = 1;
direction = 1;


end

