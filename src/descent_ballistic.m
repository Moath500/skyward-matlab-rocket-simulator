function [dY, parout] = descent_ballistic(t, Y, settings, uw, vw, ww, uncert, Hour, Day)
%{ 

ASCENT - ode function of the 6DOF Rigid Rocket Model

INPUTS:      
            - t, integration time;
            - Y, state vector, [ x y z | u v w | p q r | q0 q1 q2 q3 ]:

                                * (x y z), NED{north, east, down} horizontal frame; 
                                * (u v w), body frame velocities;
                                * (p q r), body frame angular rates;
                                * (q0 q1 q2 q3), attitude unit quaternion.

            - settings, rocket data structure;
            - uw, wind component along x;
            - vw, wind component along y;
            - ww, wind component along z;
            - uncert, wind uncertanties;
            - Hour, hour of the day of the needed simulation;
            - Day, day of the month of the needed simulation;

OUTPUTS:    
            - dY, state derivatives;
            - parout, interesting fligth quantities structure (aerodyn coefficients, forces and so on..).


NOTE: To get the NED velocities the body-frame must be multiplied for the
conjugated of the current attitude quaternion
E.G.  quatrotate(quatconj(Y(:,10:13)),Y(:,4:6))


Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 31.XII.2014

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 13/01/2018

%}

% recalling the state
x = Y(1);
y = Y(2);
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);
p = Y(7);
q = Y(8);
r = Y(9);
q0 = Y(10);
q1 = Y(11);
q2 = Y(12);
q3 = Y(13);
m = settings.m0;
Ixx = settings.Ixxe;
Iyy = settings.Iyye;
Izz = settings.Izze;

[lat, lon, ~] = ned2geodetic(x, y, 0, settings.lat0, settings.lon0, 0, wgs84Ellipsoid);     % geographic coordinates


Q = [ q0 q1 q2 q3];
Q_conj = [ q0 -q1 -q2 -q3];
normQ = norm(Q);

Q = Q/normQ;

%% ADDING WIND (supposed to be added in NED axes);

if settings.wind.model
   
    if settings.stoch.N > 1
        [uw,vw,ww] = wind_matlab_generator(settings,z,t,Hour,Day);
    else
        [uw,vw,ww] = wind_matlab_generator(settings,z,t);
    end
    
elseif settings.wind.input
    
    [uw,vw,ww] = wind_input_generator(settings,z,uncert);
    
end

    wind = quatrotate(Q, [uw vw ww]);

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
Vels = quatrotate(Q_conj,[u v w]);
V_norm = norm([ur vr wr]);

%% CONSTANTS
% Everything related to empty condition (descent-fase)

S = settings.S;              % [m^2] cross surface
C = settings.C;              % [m]   caliber
CoeffsE = settings.CoeffsE;  % [/] Empty Rocket Coefficients
g = 9.80655;                 % [N/kg] module of gravitational field at zero
T = 0;
    
%% ATMOSPHERE DATA

if -z < 0     % z is directed as the gravity vector
    z = 0;
end

[~, a, P, rho] = atmosisa(-z+settings.z0);
M = V_norm/a;
M_value = M;

%% AERODYNAMICS ANGLES

if not(u < 1e-3 || V_norm < 1e-3)
    alpha = atan(wr/ur);
    beta = asin(vr/V_norm);
else
    alpha = 0;
    beta = 0;
end

alpha_value = alpha;
beta_value = beta;


%% DATCOM COEFFICIENTS

givA = settings.Alphas*pi/180;
givB = settings.Betas*pi/180;
givH = settings.Altitudes;
givM = settings.Machs;

%% INTERPOLATION AT THE BOUNDARIES

if M > givM(end)
   
    M = givM(end);

end

if M < givM(1)
   
    M = givM(1);

end

if alpha > givA(end)
    
    alpha = givA(end);

elseif alpha < givA(1)
       
        alpha = givA(1);

end

if beta > givB(end)
    beta = givB(end);
elseif beta < givB(1)
        beta = givB(1);
end

if -z > givH(end)
    z = -givH(end);
    
elseif -z < givH(1)
        z = -givH(1);
end

%% CHOSING THE CONDITION VALUE
% interpolation of the coefficients with the value in the nearest condition of the Coeffs matrix

CA = interp4_easy(givA,givM,givB,givH,CoeffsE.CA,alpha,M,beta,-z);
CYB = interp4_easy(givA,givM,givB,givH,CoeffsE.CYB,alpha,M,beta,-z);
CNA = interp4_easy(givA,givM,givB,givH,CoeffsE.CNA,alpha,M,beta,-z);
Cl = interp4_easy(givA,givM,givB,givH,CoeffsE.CLL,alpha,M,beta,-z);
Clp = interp4_easy(givA,givM,givB,givH,CoeffsE.CLLP,alpha,M,beta,-z);
Cma = interp4_easy(givA,givM,givB,givH,CoeffsE.CMA,alpha,M,beta,-z);
Cmad = interp4_easy(givA,givM,givB,givH,CoeffsE.CMAD,alpha,M,beta,-z);
Cmq = interp4_easy(givA,givM,givB,givH,CoeffsE.CMQ,alpha,M,beta,-z);
Cnb = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNB,alpha,M,beta,-z);
Cnr = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNR,alpha,M,beta,-z);
Cnp = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNP,alpha,M,beta,-z);

%% FORCES
% first computed in the body-frame reference system

qdyn = 0.5*rho*V_norm^2;        % [Pa] dynamics pressure
qdynL_V = 0.5*rho*V_norm*S*C;   % 

X = qdyn*S*CA;                  % [N] x-body component of the aerodynamics force
Y = qdyn*S*CYB*beta;            % [N] y-body component of the aerodynamics force
Z = qdyn*S*CNA*alpha;           % [N] z-body component of the aerodynamics force
Fg = quatrotate(Q,[0 0 m*g])';  % [N] force due to the gravity

F = Fg +[-X,+Y,-Z]';            % [N] total forces vector

%% STATE DERIVATIVES

% velocity
du = F(1)/m-q*w+r*v;
dv = F(2)/m-r*u+p*w;
dw = F(3)/m-p*v+q*u;

% Rotation
dp = (Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*Cl+Clp*p*C/2);
dq = (Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cma*alpha + (Cmad+Cmq)*q*C/2);
dr = (Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cnb*beta + (Cnr*r+Cnp*p)*C/2);

% Quaternion
OM = 1/2* [ 0 -p -q -r  ;
            p  0  r -q  ;
            q -r  0  p  ;
            r  q -p  0 ];

dQQ = OM*Q'; 

%% FINAL DERIVATIVE STATE ASSEMBLING

dY(1:3) = Vels;
dY(4) = du;
dY(5) = dv;
dY(6) = dw;
dY(7) = dp;
dY(8) = dq;
dY(9) = dr;
dY(10:13) = dQQ;
dY=dY';

%% SAVING THE QUANTITIES FOR THE PLOTS

parout.integration.t = t;

parout.interp.M = M_value;
parout.interp.alpha = alpha_value;
parout.interp.beta = beta_value;
parout.interp.alt = -z;

parout.wind.NED_wind = [uw, vw, ww];
parout.wind.body_wind = wind;

parout.forces.AeroDyn_Forces = [X, Y, Z];
parout.forces.T = T;

parout.air.rho = rho;
parout.air.P = P;

parout.accelerations.body_acc = [du, dv, dw];
parout.accelerations.ang_acc = [dp, dq, dr];

parout.forces.AeroDyn_Forces = [X, Y, Z];
parout.forces.T = T;
parout.coeff.CA = CA;
parout.coeff.CYB = CYB;
parout.coeff.CNA = CNA;
parout.coeff.Cl = Cl;
parout.coeff.Clp = Clp;
parout.coeff.Cma = Cma;
parout.coeff.Cmad = Cmad;
parout.coeff.Cmq = Cmq;
parout.coeff.Cnb = Cnb;
parout.coeff.Cnr = Cnr;
parout.coeff.Cnp = Cnp;

parout.geo_cord.lat = lat;
parout.geo_cord.lon = lon;

end