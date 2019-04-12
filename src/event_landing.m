function [value,isterminal,direction] = event_landing(~,Y,settings,varargin)
% Event Function for ODE to stop simulation @ landing

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

x = Y(1);
y = Y(2);
z = -Y(3);

if settings.rocket_name == "R2A_hermes" && settings.terrain
    zloc = -settings.funZ(x,y);
    if zloc > 853
        zloc = 853;
    end
    
    if zloc < -656
        zloc = -656;
    end
    
    value = z - zloc;
else
    value = z;
end



isterminal = 1;
direction = 0;


end

