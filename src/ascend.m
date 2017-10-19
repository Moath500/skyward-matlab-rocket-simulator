function [ dY ] = ascend( t,Y,settings,uw,vw,ww )
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

global bool
% if bool = 0 the trends during the integration are NOT requested
% if bool = 1 the trends during the integration are saved and plotted
global t_plot contatore
if bool == 1
    t_plot(contatore) = t;
end


Q = [ q0 q1 q2 q3];
Q_conj= [ q0 -q1 -q2 -q3];
normQ=norm(Q);

if abs(normQ-1)>0.1
    Q = Q/normQ;
    %frprintf('Norm of the quaternion is (too much!) less than one (%3.2f)\n', ...
        %normQ);
end


% Adding Wind (supposed to be added in NED axes);

if settings.wind.model
% Wind Model
[uw,vw,ww] = wind_generator(settings,z,t,Q);
wind = [ uw,vw,ww ];
else
% constant wind
wind = quatrotate(Q, [uw vw ww]);
end

% % Wind model
% h = 0;
% if -z > 0
%     h = -z+settings.z0;
% end
% % global alt
% % if bool == 1
% %     alt = [alt; h];
% % end
%
% % wind in NED axis
% [uw,vw] = atmoshwm07(settings.wind.Lat, settings.wind.Long, h, ...
%     'day',settings.wind.Day,'model','quiet');
%
% % global WIND
% % if bool == 1
% %     WIND = [WIND; uw, vw, ww];
% % end
%
% % wind in body frame
% wind = quatrotate(Q, [uw, vw, ww]);

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
Vels = quatrotate(Q_conj,[u v w]);

V_norm = norm([ur vr wr]);



% Constants
S = settings.S;
C = settings.C;
CoeffsE = settings.CoeffsE; %Empty Rocket Coefficients
CoeffsF = settings.CoeffsF; % Full Rocket Coefficients

g = 9.81;
%Time-Dependant Properties
m0 = settings.m0; %kg mass
mfr = settings.mfr;
tb = settings.tb;

Ixxf= settings.Ixxf;
Ixxe= settings.Ixxe;
Iyyf= settings.Iyyf;
Iyye= settings.Iyye;
Izzf= settings.Izzf;
Izze= settings.Izze;

dI = 1/tb*([Ixxf Iyyf Izzf]'-[Ixxe Iyye Izze]');

if t<tb
    mdot = -mfr;
    Ixxdot = -dI(1);
    Iyydot = -dI(2);
    Izzdot = -dI(3);
%     T = settings.T_coeffs*[t^0 t^1 t^2 t^3 t^4]';
    T = interp1(settings.motor.exp_time, settings.motor.exp_thrust, t);

else
    mdot = 0;
    Ixxdot = 0;
    Iyydot = 0;
    Izzdot = 0;
    T=0;
end

global T_plot
if bool == 1
    T_plot(contatore) = T; % thrust
end

%Atmosphere
if -z<0
    z = 0;
end
[~, a, ~, rho] = atmoscoesa(-z+settings.z0);
M = V_norm/a;

global M_plot
if bool == 1
    M_plot(contatore) = M; % mach
end

%Aero Angles
if not(u<1e-1 || V_norm<1e-3);
    alpha = atan(wr/ur);
    beta = asin(vr/V_norm);
else
    alpha =0;
    beta = 0;
end

global alpha_plot beta_plot
if bool == 1
    alpha_plot(contatore) = alpha;
    beta_plot(contatore) = beta;
end

% beta = 0; % prova simulazione con beta imposto uguale a 0



% Coeffs. Interpolation
givA = settings.Alphas*pi/180;
givB = settings.Betas*pi/180;
givH = settings.Altitudes;
givM = settings.Machs;

% Interpolation error check

if M > givM(end)
    %frprintf('Mach No. is more than the max provided (M = %3.2f @ t=%2.2f)\n\n', ...
        %M,t);
    M = givM(end);

end

if M<givM(1)
    %frprintf('Mach No. is less than the min provided (M = %3.2f @ t=%2.2f)\n\n', ...
       % M,t);
    M = givM(1);

end

if alpha > givA(end)
    %frprintf('AoA is more than the max provided (alpha = %3.2f @ t=%2.2f)\n\n', ...
            %alpha*180/pi,t)
    alpha= givA(end);

else if alpha < givA(1)
        %frprintf('AoA is less than the min provided (alpha = %3.2f @ t=%2.2f)\n\n',...
            %alpha*180/pi,t);
        alpha = givA(1);


    end
end

if beta > givB(end)
    beta= givB(end);
else if beta < givB(1)
        beta = givB(1);
    end
end

if -z >givH(end)
    z = -givH(end);
else if -z < givH(1)
        z = -givH(1);
    end
end
global XCP
% interpolation with the value in the nearest condition: interp_easy
% full
CAf=interp4_easy(givA,givM,givB,givH,CoeffsF.CA,alpha,M,beta,-z);%,'nearest');
CYBf=interp4_easy(givA,givM,givB,givH,CoeffsF.CYB,alpha,M,beta,-z);%,'nearest');
CNAf=interp4_easy(givA,givM,givB,givH,CoeffsF.CNA,alpha,M,beta,-z);%,'nearest');
Clf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLL,alpha,M,beta,-z);%,'nearest');
Clpf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLLP,alpha,M,beta,-z);%,'nearest');
Cmaf=interp4_easy(givA,givM,givB,givH,CoeffsF.CMA,alpha,M,beta,-z);%,'nearest');
Cmadf=interp4_easy(givA,givM,givB,givH,CoeffsF.CMAD,alpha,M,beta,-z);%,'nearest');
Cmqf=interp4_easy(givA,givM,givB,givH,CoeffsF.CMQ,alpha,M,beta,-z);%,'nearest');
Cnbf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLNB,alpha,M,beta,-z);%,'nearest');
Cnrf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLNR,alpha,M,beta,-z);%,'nearest');
Cnpf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLNP,alpha,M,beta,-z);%,'nearest');
%
Clrf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLLR,alpha,M,beta,-z);%,'nearest');
Clbf=interp4_easy(givA,givM,givB,givH,CoeffsF.CLLB,alpha,M,beta,-z);%,'nearest');
CDf =interp4_easy(givA,givM,givB,givH,CoeffsF.CD,alpha,M,beta,-z);
XCPf =interp4_easy(givA,givM,givB,givH,CoeffsF.X_C_P,alpha,M,beta,-z);
% empty
CAe=interp4_easy(givA,givM,givB,givH,CoeffsE.CA,alpha,M,beta,-z);%,'nearest');
CYBe=interp4_easy(givA,givM,givB,givH,CoeffsE.CYB,alpha,M,beta,-z);%,'nearest');
CNAe=interp4_easy(givA,givM,givB,givH,CoeffsE.CNA,alpha,M,beta,-z);%,'nearest');
Cle=interp4_easy(givA,givM,givB,givH,CoeffsE.CLL,alpha,M,beta,-z);%,'nearest');
Clpe=interp4_easy(givA,givM,givB,givH,CoeffsE.CLLP,alpha,M,beta,-z);%,'nearest');
Cmae=interp4_easy(givA,givM,givB,givH,CoeffsE.CMA,alpha,M,beta,-z);%,'nearest');
Cmade=interp4_easy(givA,givM,givB,givH,CoeffsE.CMAD,alpha,M,beta,-z);%,'nearest');
Cmqe=interp4_easy(givA,givM,givB,givH,CoeffsE.CMQ,alpha,M,beta,-z);%,'nearest');
Cnbe=interp4_easy(givA,givM,givB,givH,CoeffsE.CLNB,alpha,M,beta,-z);%,'nearest');
Cnre=interp4_easy(givA,givM,givB,givH,CoeffsE.CLNR,alpha,M,beta,-z);%,'nearest');
Cnpe=interp4_easy(givA,givM,givB,givH,CoeffsE.CLNP,alpha,M,beta,-z);%,'nearest');
%
Clre=interp4_easy(givA,givM,givB,givH,CoeffsE.CLLR,alpha,M,beta,-z);%,'nearest');
Clbe=interp4_easy(givA,givM,givB,givH,CoeffsE.CLLB,alpha,M,beta,-z);%,'nearest');
CDe =interp4_easy(givA,givM,givB,givH,CoeffsE.CD,alpha,M,beta,-z);
XCPe =interp4_easy(givA,givM,givB,givH,CoeffsE.X_C_P,alpha,M,beta,-z);
% % linear interpolation of the coeff
% %full
% CAf = interpn(givA,givM,givB,givH,CoeffsF.CA,alpha,M,beta,-z);
% CYBf = interpn(givA,givM,givB,givH,CoeffsF.CYB,alpha,M,beta,-z);
% CNAf = interpn(givA,givM,givB,givH,CoeffsF.CNA,alpha,M,beta,-z);
% Clf = interpn(givA,givM,givB,givH,CoeffsF.CLL,alpha,M,beta,-z);
% Clpf = interpn(givA,givM,givB,givH,CoeffsF.CLLP,alpha,M,beta,-z);
% Cmaf = interpn(givA,givM,givB,givH,CoeffsF.CMA,alpha,M,beta,-z);
% Cmadf = interpn(givA,givM,givB,givH,CoeffsF.CMAD,alpha,M,beta,-z);
% Cmqf = interpn(givA,givM,givB,givH,CoeffsF.CMQ,alpha,M,beta,-z);
% Cnbf = interpn(givA,givM,givB,givH,CoeffsF.CLNB,alpha,M,beta,-z);
% Cnrf = interpn(givA,givM,givB,givH,CoeffsF.CLNR,alpha,M,beta,-z);
% Cnpf = interpn(givA,givM,givB,givH,CoeffsF.CLNP,alpha,M,beta,-z);
% %
% Clrf = interpn(givA,givM,givB,givH,CoeffsF.CLLR,alpha,M,beta,-z);
% Clbf = interpn(givA,givM,givB,givH,CoeffsF.CLLB,alpha,M,beta,-z);
% CDf = interpn(givA,givM,givB,givH,CoeffsF.CD,alpha,M,beta,-z);
%
% % empty
% CAe = interpn(givA,givM,givB,givH,CoeffsE.CA,alpha,M,beta,-z);
% CYBe = interpn(givA,givM,givB,givH,CoeffsE.CYB,alpha,M,beta,-z);
% CNAe = interpn(givA,givM,givB,givH,CoeffsE.CNA,alpha,M,beta,-z);
% Cle = interpn(givA,givM,givB,givH,CoeffsE.CLL,alpha,M,beta,-z);
% Clpe = interpn(givA,givM,givB,givH,CoeffsE.CLLP,alpha,M,beta,-z);
% Cmae = interpn(givA,givM,givB,givH,CoeffsE.CMA,alpha,M,beta,-z);
% Cmade = interpn(givA,givM,givB,givH,CoeffsE.CMAD,alpha,M,beta,-z);
% Cmqe = interpn(givA,givM,givB,givH,CoeffsE.CMQ,alpha,M,beta,-z);
% Cnbe = interpn(givA,givM,givB,givH,CoeffsE.CLNB,alpha,M,beta,-z);
% Cnre = interpn(givA,givM,givB,givH,CoeffsE.CLNR,alpha,M,beta,-z);
% Cnpe = interpn(givA,givM,givB,givH,CoeffsE.CLNP,alpha,M,beta,-z);
% %
% Clre = interpn(givA,givM,givB,givH,CoeffsE.CLLR,alpha,M,beta,-z);
% Clbe = interpn(givA,givM,givB,givH,CoeffsE.CLLB,alpha,M,beta,-z);
% CDe =interpn(givA,givM,givB,givH,CoeffsE.CD,alpha,M,beta,-z);

%Linear interpolation from empty and full configuration coefficients
if t<tb
    CA = t/tb*(CAe-CAf)+CAf;
    CYB= t/tb*(CYBe-CYBf)+CYBf;
    CNA= t/tb*(CNAe-CNAf)+CNAf;
    Cl = t/tb*(Cle-Clf)+Clf;
    Clp= t/tb*(Clpe-Clpf)+Clpf;
    Cma= t/tb*(Cmae-Cmaf)+Cmaf;
    Cmad=t/tb*(Cmade-Cmadf)+Cmadf;
    Cmq= t/tb*(Cmqe-Cmqf)+Cmqf;
    Cnb= t/tb*(Cnbe-Cnbf)+Cnbf;
    Cnr= t/tb*(Cnre-Cnrf)+Cnrf;
    Cnp= t/tb*(Cnpe-Cnpf)+Cnpf;
    %
    Clr= t/tb*(Clre-Clrf)+Clrf;
    Clb= t/tb*(Clbe-Clbf)+Clbf;
    CD= t/tb*(CDe-CDf)+CDf;
    XCP(contatore)= t/tb*(XCPe-XCPf)+XCPf;
else
    CA = CAe;
    CYB= CYBe;
    CNA= CNAe;
    Cl = Cle;
    Clp= Clpe;
    Cma= Cmae;
    Cmad=Cmade;
    Cmq= Cmqe;
    Cnb= Cnbe;
    Cnr= Cnre;
    Cnp= Cnpe;
    %
    Clr= Clre;
    Clb= Clbe;
    CD = CDe;
    XCP(contatore) = XCPe;
end


global CA_plot % coefficiente aerodinamico
if bool == 1
    CA_plot(contatore) = CD; % Axial force coeff
end


%Center of Mass



%Body Frame

qdynS=0.5*rho*V_norm^2*S;
qdynL_V = 0.5*rho*V_norm*S*C;

X=qdynS*CA;
Y=qdynS*CYB*beta;
Z=qdynS*CNA*alpha;

global Drag_plot
if bool == 1
    Drag_plot(contatore) = qdynS*CD;
end


%Forces
F = quatrotate(Q,[0 0 m*g])';

F = F +[-X+T,+Y,-Z]';

global Forces_plot
if bool == 1
    Forces_plot(contatore) = F(1);
end


du=F(1)/m -q*w + r*v;% - mdot/m*u;
dv=F(2)/m -r*u + p*w;% - mdot/m*v;
dw=F(3)/m -p*v + q*u;% - mdot/m*w;




% Rotation
dp=(Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*Cl+Clp*p*C/2)-Ixxdot*p/Ixx;
dq=(Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cma*alpha + (Cmad+Cmq)*q*C/2)...
    -Iyydot*q/Iyy;
dr=(Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cnb*beta + (Cnr*r+Cnp*p)*C/2)...
    -Izzdot*r/Izz;


% dp=(Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*(Cl+Clb*beta)+(Clp*p+Clr*r)*C/2)...
%     -Ixxdot*p/Ixx;
% dq=(Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cma*alpha + (Cmad+Cmq)*q*C/2)...
%     -Iyydot*q/Iyy;
% dr=(Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cnb*beta + (Cnr*r+Cnp*p)*C/2)...
%     -Izzdot*r/Izz;



% Launch Pad-relative coordinates
Xb = quatrotate(Q,[x y z]);

%

%On Launch Pad No Torque
 if Xb(1) < settings.lrampa
    dp = 0;
    dq = 0;
    dr = 0;

    %Ground Reaction.
    if T < m*g
        du = 0;
        dv = 0;
        dw = 0;

    end
end







OM = 1/2*[
   0 -p -q -r
   p 0 r -q
   q -r 0 p
   r q -p 0];


dQQ = OM*Q';

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
dY=dY';

if bool == 1
    contatore = contatore + 1;
end

end
