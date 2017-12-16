function [dY] = descent_parachute(t,Y,settings,uw,vw,ww,para)
% ODE-Function for Parachute descent 
% State = ( x y z | u v w  )

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

% x = Y(1);
% y = Y(2);
  z = Y(3);
  u = Y(4);
  v = Y(5);
  w = Y(6);

%% ADDING WIND (supposed to be added in NED axes);

if settings.wind.model
    [uw,vw,ww] = wind_matlab_generator(settings,z,t);
    wind = [uw,vw,ww];

else

wind = [uw vw ww]; % constant wind

end

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

V_norm = norm([ur vr wr]);

%% CONSTANTS

switch para
    case 1
        S = settings.para1.S;                                               % [m^2]   Surface
        CD = settings.para1.CD;                                             % [/] Parachute Drag Coefficient
        CL = settings.para1.CL;                                             % [/] Parachute Lift Coefficient
        pmass = 0;                                                          % [kg] detached mass
    case 2
        S = settings.para2.S;                                               % [m^2]   Surface
        CD = settings.para2.CD;                                             % [/] Parachute Drag Coefficient
        CL = settings.para2.CL;                                             % [/] Parachute Lift Coefficient
        pmass = settings.para1.mass + settings.mnc;                         % [kg] detached mass(drogue1 + nosecone)
    case 3
        S = settings.para3.S;                                               % [m^2]   Surface
        CD = settings.para3.CD;                                             % [/] Parachute Drag Coefficient
        CL = settings.para3.CL;                                             % [/] Parachute Lift Coefficient
        pmass = settings.para1.mass + settings.para2.mass + settings.mnc;   % [kg] detached mass(drogue1/2 + nosecone)
    otherwise
end

g = 9.80655;                                                                % [N/kg] magnitude of the gravitational field at zero
m = settings.ms - pmass;                                                    % [kg] descend mass

%% ATMOSPHERE DATA

if -z < 0
    z = 0;
end

[~, ~, ~, rho] = atmoscoesa(-z);


%% REFERENCE FRAME
% The parachutes are approximated as rectangular surfaces with the normal
% vector perpendicular to the relative velocity

t_vect = [ur vr wr];                     % Tangenzial vector
h_vect = [-vr ur 0];                     % horizontal vector

if (h_vect(1) == 0) && (h_vect(2) == 0)
    h_vect = [-v u 0];
end

t_vers = t_vect/norm(t_vect);            % Tangenzial versor
h_vers = -h_vect/norm(h_vect);           % horizontal versor
n_vect = cross(t_vers, h_vers);          % Normal vector
n_vers = n_vect/norm(n_vect);            % Normal versor

if (n_vers(3) > 0)                       % If the normal vector is downward directed
    n_vect = cross(h_vers, t_vers); 
    n_vers = n_vect/norm(n_vect);
end

%% FORCES

D = 0.5*rho*V_norm^2*S*CD*t_vers';       % [N] Drag vector
L = 0.5*rho*V_norm^2*S*CL*n_vers';       % [N] Lift vector
Fg = m*g*[0 0 1]';                       % [N] Gravitational Force vector                  
F = -D+L+Fg;                             % [N] total forces vector

%% STATE DERIVATIVES

% velocity
du = F(1)/m;
dv = F(2)/m;
dw = F(3)/m;

%% FINAL DERIVATIVE STATE ASSEMBLING

dY(1:3) = [u v w]';
dY(4) = du;
dY(5) = dv;
dY(6) = dw;

dY = dY';

end