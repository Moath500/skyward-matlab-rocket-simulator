function [ value,isterminal,direction ] = main_opening( t,Y,settings,varargin)
%Event that sets the main parachute opening altitude

value = Y(3)+settings.zmain;
isterminal = 1;
direction = 1;


end

