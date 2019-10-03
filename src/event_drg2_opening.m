function [value,isterminal,direction] = event_drg2_opening(~,Y,settings,varargin)
% Event that sets the main parachute opening altitude

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

x = Y(1);
y = Y(2);
z = -Y(3);

% if settings.rocket_name == "R2A_hermes" && settings.terrain
if settings.terrain
    zloc = -settings.funZ(x,y);
    if zloc > 859
        zloc = 859;
    end
    
    if zloc < -845
        zloc = -845;
    end
    
    value = z - zloc - settings.zdrg2;
else
    value = z - settings.zdrg2;
end

isterminal = 1;
direction = 0;

end

