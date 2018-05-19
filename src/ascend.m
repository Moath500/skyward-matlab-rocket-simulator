function [dY] = ascend(t,Y,settings,uw,vw,ww)
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
m = Y(14);
Ixx = Y(15);
Iyy = Y(16);
Izz = Y(17);


%% QUATERION ATTITUDE

Q = [q0 q1 q2 q3];
Q_conj = [q0 -q1 -q2 -q3];
normQ = norm(Q);

if abs(normQ-1) > 0.1
    Q = Q/normQ;
end


%% ADDING WIND (supposed to be added in NED axes);

if settings.wind.model
    [uw,vw,ww] = wind_matlab_generator(settings,z,t,Q);
    wind = [uw,vw,ww];
    
else
    
    wind = quatrotate(Q, [uw vw ww]); % constant wind
    
end

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
Vels = quatrotate(Q_conj,[u v w]);
V_norm = norm([ur vr wr]);

%% ATMOSPHERE DATA

if -z < 0     % z is directed as the gravity vector
    z = 0;
end
[~, a, ~, rho] = atmosisa(-z+settings.z0);
M = V_norm/a;
M_value = M;

%% CONSTANTS

S = settings.S;              % [m^2] cross surface
C = settings.C;              % [m]   caliber
CoeffsE = settings.CoeffsE;  % Empty Rocket Coefficients
CoeffsF = settings.CoeffsF;  % Full Rocket Coefficients
g = 9.80655;                 % [N/kg] module of gravitational field at zero
tb = settings.tb;            % [s]     Burning Time
mfr = settings.mfr;          % [kg/s]  Mass Flow Rate
OMEGA = settings.OMEGA;      %[rad] Elevation Angle in the launch pad

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
Ixxf = settings.Ixxf;        % [kg*m^2] Inertia to x-axis
Iyyf = settings.Iyyf;        % [kg*m^2] Inertia to y-axis
Izzf = settings.Izzf;        % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
Ixxe = settings.Ixxe;        % [kg*m^2] Inertia to x-axis
Iyye = settings.Iyye;        % [kg*m^2] Inertia to y-axis
Izze = settings.Izze;        % [kg*m^2] Inertia to z-axis


%% TIME-DEPENDENTS VARIABLES

dI = 1/tb*([Ixxf Iyyf Izzf]'-[Ixxe Iyye Izze]');

if t<tb
    mdot = -mfr;
    Ixxdot = -dI(1);
    Iyydot = -dI(2);
    Izzdot = -dI(3);
    T = interp1(settings.motor.exp_time, settings.motor.exp_thrust, t);
    
else             % for t >= tb the fligth condition is the empty one(no interpolation needed)
    mdot = 0;
    Ixxdot = 0;
    Iyydot = 0;
    Izzdot = 0;
    T=0;
end

T_value = T;

%% AERODYNAMICS ANGLES

if not(u < 1e-1 || V_norm < 1e-3)
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

%% CHOSING THE FULL CONDITION VALUE
% interpolation of the coefficients with the value in the nearest condition of the Coeffs matrix

CAf = interp4_easy(givA,givM,givB,givH,CoeffsF.CA,alpha,M,beta,-z);
CYBf = interp4_easy(givA,givM,givB,givH,CoeffsF.CYB,alpha,M,beta,-z);
CNAf = interp4_easy(givA,givM,givB,givH,CoeffsF.CNA,alpha,M,beta,-z);
Clf = interp4_easy(givA,givM,givB,givH,CoeffsF.CLL,alpha,M,beta,-z);
Clpf = interp4_easy(givA,givM,givB,givH,CoeffsF.CLLP,alpha,M,beta,-z);
Cmaf = interp4_easy(givA,givM,givB,givH,CoeffsF.CMA,alpha,M,beta,-z);
Cmadf = interp4_easy(givA,givM,givB,givH,CoeffsF.CMAD,alpha,M,beta,-z);
Cmqf = interp4_easy(givA,givM,givB,givH,CoeffsF.CMQ,alpha,M,beta,-z);
Cnbf = interp4_easy(givA,givM,givB,givH,CoeffsF.CLNB,alpha,M,beta,-z);
Cnrf = interp4_easy(givA,givM,givB,givH,CoeffsF.CLNR,alpha,M,beta,-z);
Cnpf = interp4_easy(givA,givM,givB,givH,CoeffsF.CLNP,alpha,M,beta,-z);
CDf = interp4_easy(givA,givM,givB,givH,CoeffsF.CD,alpha,M,beta,-z);
XCPf = interp4_easy(givA,givM,givB,givH,CoeffsF.X_C_P,alpha,M,beta,-z);

%% CHOSING THE EMPTY CONDITION VALUE
% interpolation of the coefficients with the value in the nearest condition of the Coeffs matrix

CAe = interp4_easy(givA,givM,givB,givH,CoeffsE.CA,alpha,M,beta,-z);
CYBe = interp4_easy(givA,givM,givB,givH,CoeffsE.CYB,alpha,M,beta,-z);
CNAe = interp4_easy(givA,givM,givB,givH,CoeffsE.CNA,alpha,M,beta,-z);
Cle = interp4_easy(givA,givM,givB,givH,CoeffsE.CLL,alpha,M,beta,-z);
Clpe = interp4_easy(givA,givM,givB,givH,CoeffsE.CLLP,alpha,M,beta,-z);
Cmae = interp4_easy(givA,givM,givB,givH,CoeffsE.CMA,alpha,M,beta,-z);
Cmade = interp4_easy(givA,givM,givB,givH,CoeffsE.CMAD,alpha,M,beta,-z);
Cmqe = interp4_easy(givA,givM,givB,givH,CoeffsE.CMQ,alpha,M,beta,-z);
Cnbe = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNB,alpha,M,beta,-z);
Cnre = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNR,alpha,M,beta,-z);
Cnpe = interp4_easy(givA,givM,givB,givH,CoeffsE.CLNP,alpha,M,beta,-z);
CDe = interp4_easy(givA,givM,givB,givH,CoeffsE.CD,alpha,M,beta,-z);
XCPe = interp4_easy(givA,givM,givB,givH,CoeffsE.X_C_P,alpha,M,beta,-z);

%% LINEAR INTERPOLATION BETWEEN THE TWO CONDITIONS
% Computing the value of the aerodynamics coefficients at a certain time
% Needed only for t<tb because for t>=tb the condition is the empty one

if t < tb
    CA = t/tb*(CAe-CAf)+CAf;
    CYB = t/tb*(CYBe-CYBf)+CYBf;
    CNA = t/tb*(CNAe-CNAf)+CNAf;
    Cl = t/tb*(Cle-Clf)+Clf;
    Clp = t/tb*(Clpe-Clpf)+Clpf;
    Cma = t/tb*(Cmae-Cmaf)+Cmaf;
    Cmad = t/tb*(Cmade-Cmadf)+Cmadf;
    Cmq = t/tb*(Cmqe-Cmqf)+Cmqf;
    Cnb = t/tb*(Cnbe-Cnbf)+Cnbf;
    Cnr = t/tb*(Cnre-Cnrf)+Cnrf;
    Cnp = t/tb*(Cnpe-Cnpf)+Cnpf;
    CD = t/tb*(CDe-CDf)+CDf;
    XCP_value = t/tb*(XCPe-XCPf)+XCPf;
else
    CA = CAe;
    CYB = CYBe;
    CNA = CNAe;
    Cl = Cle;
    Clp = Clpe;
    Cma = Cmae;
    Cmad =Cmade;
    Cmq = Cmqe;
    Cnb = Cnbe;
    Cnr = Cnre;
    Cnp = Cnpe;
    CD = CDe;
    XCP_value = XCPe;
end

if -z < settings.lrampa*sin(OMEGA)      % No torque on the Launch
    
    Fg = m*g*sin(OMEGA);                %[N] force due to the gravity
    X = 0.5*rho*V_norm^2*S*CA;
    F = -Fg +T -X;
    du = F/m;
    dv = 0;
    dw = 0;
    dp = 0;
    dq = 0;
    dr = 0;
    alpha_value = NaN;
    beta_value = 0;
    qdyn = NaN;
    Y = NaN;
    Z = NaN;
    XCP_value = 0;
    
    
    if T < Fg                           % No velocity untill T = Fg
        du = 0;
    end
    
else
    
    %% FORCES
    % first computed in the body-frame reference system
    
    qdyn = 0.5*rho*V_norm^2;        %[Pa] dynamics pressure
    qdynL_V = 0.5*rho*V_norm*S*C;   %
    
    X = qdyn*S*CA;                  %[N] x-body component of the aerodynamics force
    Y = qdyn*S*CYB*beta;            %[N] y-body component of the aerodynamics force
    Z = qdyn*S*CNA*alpha;           %[N] z-body component of the aerodynamics force
    Fg = quatrotate(Q,[0 0 m*g])';  %[N] force due to the gravity
    
    F = Fg +[-X+T,+Y,-Z]';          %[N] total forces vector
    
    %% STATE DERIVATIVES
    
    % velocity
    du = F(1)/m-q*w+r*v;
    dv = F(2)/m-r*u+p*w;
    dw = F(3)/m-p*v+q*u;
    
    % Rotation
    dp = (Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*Cl+Clp*p*C/2)-Ixxdot*p/Ixx;
    dq = (Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cma*alpha + (Cmad+Cmq)*q*C/2)...
        -Iyydot*q/Iyy;
    dr = (Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cnb*beta + (Cnr*r+Cnp*p)*C/2)...
        -Izzdot*r/Izz;
    
end
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
dY(14) = mdot;
dY(15) = Ixxdot;
dY(16) = Iyydot;
dY(17) = Izzdot;
dY = dY';



%% PERSISTENT VARIABLES

persistent XCP contatore t_plot ...
    T_plot M_plot CA_plot alpha_plot beta_plot Drag_plot Forces_plot F_aero_x F_aero_y F_aero_z


%% SAVING THE QUANTITIES FOR THE PLOTS

if settings.plots
    
    if isempty (contatore)
        contatore = 1;
        t_plot(contatore) = 0;
        T_plot(contatore) = 0;
        CA_plot(contatore) = 0;
        Drag_plot(contatore) = 0;
        Forces_plot(contatore) = 0;
        M_plot(contatore) = 0;
        XCP(contatore) = 0;
        alpha_plot(contatore) = 0;
        beta_plot(contatore) = 0;
        F_aero_x(contatore) = 0;
        F_aero_y(contatore) = 0;
        F_aero_z(contatore) = 0;
        
    end
    t_plot(contatore) = t;
    beta_plot(contatore) = beta_value;
    alpha_plot(contatore) = alpha_value;
    XCP(contatore) = XCP_value;
    M_plot(contatore) = M_value;
    T_plot(contatore) = T_value;
    CA_plot(contatore) = CD;
    Drag_plot(contatore) = qdyn*S*CD;
    Forces_plot(contatore) = F(1);
    F_aero_x(contatore) = X;
    F_aero_y(contatore) = Y;
    F_aero_z(contatore) = Z;
    contatore = contatore + 1;
    
    
    if Vels(3) <= 0 && t > tb
        
        ascend.XCP = XCP;
        ascend.t = t_plot;
        ascend.T = T_plot;
        ascend.M = M_plot;
        ascend.CA = CA_plot;
        ascend.alpha = alpha_plot;
        ascend.beta = beta_plot;
        ascend.Drag = Drag_plot;
        ascend.Forces = Forces_plot;
        ascend.Dyn_Forces(1,:) = F_aero_x; 
        ascend.Dyn_Forces(2,:) = F_aero_y;
        ascend.Dyn_Forces(3,:) = F_aero_z;
        
        if settings.stoch.N == 1
            save ('ascend_plot.mat', 'ascend')
        end
    end
end

