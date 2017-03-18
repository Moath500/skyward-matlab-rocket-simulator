function [ dY ] = descent_parachute( t,Y,settings,uw,vw,ww,para )
% ODEFUN for Parachute descent
% State = ( x y z | u v w  )

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

x = Y(1);
y = Y(2);
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);

%Adding wind

ur = u - uw;
vr = v - vw;
wr = w - ww;

V_norm = norm([ur vr wr]);

% Constants
if para == 1
    S = settings.para1.S; %Parachute Surface
    CD = settings.para1.CD; %Parachute CD
    pmass = settings.para1.mass;
else
    S = settings.para2.S;
    CD = settings.para2.CD;
    pmass = settings.para2.mass;
end

rad = sqrt(S/pi); %Equivalent parachute radius;

%Side Surface 
Ss = pi*(rad/2)^2;

g = 9.81;

m = settings.ms + pmass;

%Atmosphere
if -z<0
    z = 0;
end

[~, ~, ~, rho] = atmoscoesa(-z);

%Center of Mass



%Body Frame

D=0.5*rho*V_norm^2*S*CD;

X=D*ur/V_norm;
Y=D*vr/V_norm;
Z=D*wr/V_norm -m*g;

%Forces
F = [-X,-Y,-Z]';



du=F(1)/m;
dv=F(2)/m;
dw=F(3)/m;





dY(1:3) = [u v w]';
dY(4) = du;
dY(5) = dv;
dY(6) = dw;

dY = dY';

end
