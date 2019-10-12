function [value, isterminal, direction] = event_para_cut(~, Y, settings, varargin)
%{

EVENT_PARA_CUT - Event function to stop simulation at the chosen altitude to cut the 
                 parachute, checking when a value tends to zero;
                 the value taken is to account is the vertical position respect to the chosen altitude,
                 z - zcut = 0 --> parachute cut

INPUTS:     
            - t, integration time;
            - Y, state vector, [ x y z | u v w | p q r | q0 q1 q2 q3 | m | Ixx Iyy Izz ]:

                                * (x y z), NED{north, east, down} horizontal frame; 
                                * (u v w), body frame velocities;
                                * (p q r), body frame angular rates;
                                * m , total mass;
                                * (Ixx Iyy Izz), Inertias;
                                * (q0 q1 q2 q3), attitude unit quaternion.

            - settings, rocket data structure.

OUTPUTS:        
            - isterminal, logical input to stop the integration;
            - direction, to select the sign that the function must have when tends to zero, 1 = positive;
            - value, selected value to check if the integration has to be stopped (vertical position).

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 25.IV.2014

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Last Revision: 09/10/2019

%}

x = Y(1);
y = Y(2);
z = -Y(3);

para = varargin{4};

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

