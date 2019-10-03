function [value,isterminal,direction] = event_rog_opening(~,Y,settings,varargin)
%Event that sets the main parachute opening altitude


% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

value = -Y(3)-settings.zrog;
isterminal = 1;
direction = 1;


end

