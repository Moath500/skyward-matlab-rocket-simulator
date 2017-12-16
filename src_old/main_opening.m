function [ value,isterminal,direction ] = main_opening( t,Y,settings,varargin)
%Event that sets the main parachute opening altitude


% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

value = Y(3)+settings.zmain;
isterminal = 1;
direction = 1;


end

