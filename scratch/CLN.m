close all
clear, clc

% AERODYNAMICS DETAILS %
% This coefficients are intended to be obtained through MISSILE DATCOM
% (than parsed with the tool datcom_parser.py)
CoeffsF = load('for006_full.mat','Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');

%Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
%empty and full configuration
s = load('for006_full.mat','State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
clear('s');

% Coeffs. Interpolation
givA = settings.Alphas*pi/180;
givB = settings.Betas*pi/180;
givH = settings.Altitudes;
givM = settings.Machs;

alpha = 0;
M = 1.7;
beta = [-5:0.1:5]*pi/180;
z = 0;

Cln = zeros(length(beta),2);
Cnb = zeros(length(beta),2);
Cnr = zeros(length(beta),2);
Cnp = zeros(length(beta),2);

for k = 1:length(beta)
    Cln(k,1)=interp4_easy(givA,givM,givB,givH,settings.CoeffsF.CLN,alpha,M,beta(k),-z);%,'nearest');
    Cnb(k,1)=interp4_easy(givA,givM,givB,givH,settings.CoeffsF.CLNB,alpha,M,beta(k),-z);%,'nearest');
    Cnr(k,1)=interp4_easy(givA,givM,givB,givH,settings.CoeffsF.CLNR,alpha,M,beta(k),-z);%,'nearest');
    Cnp(k,1)=interp4_easy(givA,givM,givB,givH,settings.CoeffsF.CLNP,alpha,M,beta(k),-z);%,'nearest');
    
    Cln(k,2)=interpn(givA,givM,givB,givH,settings.CoeffsF.CLN,alpha,M,beta(k),-z);%,'nearest');
    Cnb(k,2)=interpn(givA,givM,givB,givH,settings.CoeffsF.CLNB,alpha,M,beta(k),-z);%,'nearest');
    Cnr(k,2)=interpn(givA,givM,givB,givH,settings.CoeffsF.CLNR,alpha,M,beta(k),-z);%,'nearest');
    Cnp(k,2)=interpn(givA,givM,givB,givH,settings.CoeffsF.CLNP,alpha,M,beta(k),-z);%,'nearest');
end

figure(), plot(beta*180/pi, Cln), title('CLN'), grid on;
figure(), plot(beta*180/pi, Cnb), title('CLN/B'), grid on;
figure(), plot(beta*180/pi, Cnr), title('CLN/r'), grid on;
figure(), plot(beta*180/pi, Cnp), title('CLN/p'), grid on;