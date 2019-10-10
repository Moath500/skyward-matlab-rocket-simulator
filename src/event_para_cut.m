function [value, isterminal, direction] = event_para_cut(~, Y, settings, varargin)

%{
 EVENT_PARA_CUT: Event function that determines when the parachute has to be cutted

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept
email: adriano.filippo.inno@skywarder.eu
Revision date: 09/10/2019

%}

x = Y(1);
y = Y(2);
z = -Y(3);

para = varargin{4};

% if settings.rocket_name == "R2A_hermes" && settings.terrain
if settings.terrain
    zloc = -settings.funZ(x,y);
    if zloc > 859
        zloc = 859;
    end
    
    if zloc < -845
        zloc = -845;
    end
    
    value = z - zloc - settings.para(para).z_cut;
else
    value = z - settings.para(para).z_cut;
end

isterminal = 1;
direction = 0;

end

