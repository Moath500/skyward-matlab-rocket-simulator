function [ value,isterminal,direction ] = drg2_opening( t,Y,settings,varargin)
%Event that sets the main parachute opening altitude

value = Y(3)+settings.zdrg2;
isterminal = 1;
direction = 1;


end

