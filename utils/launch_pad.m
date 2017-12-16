% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

close all
clear, clc

% inclinazione rampa
theta = 85*pi/180;

% [vv]horizon = R'*[vv]corpo
% [vv]corpo = R*[vv]horizon
R = [ cos(theta), 0, sin(theta); ...
               0, 1,          0; ...
     -sin(theta), 0, cos(theta)];

% velocit� di uscita dalla rampa
Lrampa = 4.5; % m
acc = 5*9.80665; % m/s2
texit = sqrt(Lrampa/(.5*acc)); % sec
V = acc*texit % m/s
% V = 32;% m/s
vv = V*[1, 0, 0]';

% vento
% incidenza vento
dev = [-180:5:180]*pi/180;
W = sqrt(2^2+5^2) * 1.2 % m/s
% W = 10

[~, a, ~, rho] = atmoscoesa(0);

C = 0.174; % m Caliber (Fuselage Diameter)
S = 0.02378; % m2 Cross-sectional Surface
L = 4.32050; % m lunghezza missile

load for006_full.mat
% Coeffs. Interpolation
givA = State.Alphas;
givB = State.Betas;
givH = State.Altitudes;
givM = State.Machs;

norm_vr = zeros(length(dev), 1);
M = zeros(length(dev), 1);
alpha = zeros(length(dev), 1);
beta = zeros(length(dev), 1);
XCP = zeros(length(dev), 1);
CM = zeros(length(dev), 1);
Mpitch = zeros(length(dev), 1);

for n = 1:length(dev)
    ww = -W*[cos(dev(n)), sin(dev(n)), 0]'; % assi horizon
    ww = R*ww; % assi corpo

    % vel relativa al vento
    vvr = vv - ww; % assi corpo
    norm_vr(n) = norm(vvr);
    M(n) = norm_vr(n)/a;

    ur = vvr(1);
    vr = vvr(2);
    wr = vvr(3);

    % alpha
    alpha(n) = atan(-wr/ur) *180/pi; % deg
    beta(n) = asin(vr/norm(vvr)) *180/pi; % deg

    XCP(n) = interpn(givA,givM,givB,givH,Coeffs.X_C_P,alpha(n),M(n),beta(n),0);
    CM(n) = interpn(givA,givM,givB,givH,Coeffs.CM,alpha(n),M(n),beta(n),0);
    Mpitch(n) = 0.5*rho*norm_vr(n)^2*S*C*CM(n);
end

figure()
plot(dev*180/pi, XCP), grid on
xlabel('angolo vento [deg]'), ylabel('XCP [ref length]')

figure()
plot(dev*180/pi, CM), grid on
xlabel('angolo vento [deg]'), ylabel('CM [-]')

figure()
plot(dev*180/pi, Mpitch), grid on
xlabel('angolo vento [deg]'), ylabel('M pitch [Nm]')
