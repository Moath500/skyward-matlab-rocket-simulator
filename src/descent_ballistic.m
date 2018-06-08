function [dY] = descent_ballistic(t,Y,settings,uw,vw,ww,Hour,Day)
% ODE-Function of the 6DOF Rigid Rocket Model
% State = ( x y z | u v w | p q r | q0 q1 q2 q3 )
%
% (x y z): NED Earth's Surface Centered Frame ("Inertial") coordinates
% (u v w): body frame velocities
% (p q r): body frame angular rates
% (q0 q1 q2 q3): attitude unit quaternion
%
%
% NOTE: To get the NED velocities the body-frame must be multiplied for the
% conjugated of the current attitude quaternion
% E.G.
%
%
% quatrotate(quatconj(Y(:,10:13)),Y(:,4:6))

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 31.XII.2014
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


Q = [ q0 q1 q2 q3];
Q_conj = [ q0 -q1 -q2 -q3];
normQ = norm(Q);

if abs(normQ-1) > 0.1
    Q = Q/normQ;
end

%% ADDING WIND (supposed to be added in NED axes);


if settings.wind.model
    if settings.stoch.N > 1
        [uw,vw,ww] = wind_matlab_generator(settings,z,t,Hour,Day);
    else
        [uw,vw,ww] = wind_matlab_generator(settings,z,t);
    end
    wind = [uw,vw,ww];
    
else
    
    wind = [uw vw ww]; % constant wind

end

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

T = 0;                       % Thrust
    
%% ATMOSPHERE DATA

if -z < 0     % z is directed as the gravity vector
    z = 0;
end

[~, a, ~, rho] = atmoscoesa(-z+settings.z0);
M = V_norm/a;

%% AERODYNAMICS ANGLES

if not(u<1e-1 || V_norm<1e-3)
    alpha = atan(wr/ur);
    beta = asin(vr/V_norm);
else
    alpha =0;
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

F = Fg +[-X+T,+Y,-Z]';          % [N] total forces vector

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

%% PERSISTENT VARIABLES

persistent t_plot contatore beta_plot alpha_plot wind_plot alt_plot


%% SAVING THE QUANTITIES FOR THE PLOTS

if settings.plots
    
    if settings.stoch.N == 1
        
        if isempty (contatore)
            contatore = 1;
            t_plot(contatore) = 0;
            beta_plot(contatore) = 0;
            alpha_plot(contatore) = 0;
            wind_plot(:,contatore) = zeros(3,1);
            alt_plot(contatore) = 0;
        end
        
        t_plot(contatore) = t;
        beta_plot(contatore) = beta_value;
        alpha_plot(contatore) = alpha_value;
        wind_plot(:,contatore) = [uw,vw,ww];
        alt_plot(contatore) = -z;
        contatore = contatore + 1;
        
        
        descent_bal.t = t_plot;
        descent_bal.alpha = alpha_plot;
        descent_bal.beta = beta_plot;
        descent_bal.wind = wind_plot;
        descent_bal.alt = alt_plot;
        
        
        if -z <= 0
            save ('descent_plot.mat', 'descent_bal')
        end
        
        
    end
end