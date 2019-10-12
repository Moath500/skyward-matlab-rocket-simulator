function [value, isterminal, direction] = event_apogee(t, Y, settings, varargin)
%{

EVENT_APOGEE - Event function to stop simulation at apogee checking when a value tends to zero;
               the value taken is to account is the vertical velocity, vy = 0 --> apogee

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
            - value, selected value to check if the integration has to be stopped (vertical velocity).

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 25.IV.2014

%}

Q = Y(10:13)';

% Inertial Frame velocities
vels = quatrotate(quatconj(Q),Y(4:6)');

% Stop checking if I'm in Propulsion Phase
if t > settings.tb
    value = vels(3);
else
    value = 1;
end

isterminal = 1;
direction = 1;


end

