% FIN_OPT - This script evaluate a comparison between different
% aerodynamics coefficients calculated with datcom 
% 
% Author: Matteo Pozzoli
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: matteo.pozzoli@skywarder.eu
% Website: http://www.skywarder.eu
% License: 2-clause BSD
% 
% Release date: 14/10/2019

clear 
close all
clc 

%% DATA
run config.m

r_name = {'R2A' 'R1X'};
n = length(r_name);

% settings 
grafico = false;
cal_min = 1;                    % minum stability margin required

%% RUN 

% preallocation 
apogee = zeros(1,n);
max_a = zeros(1,n);
Vexit = zeros(1,n);
XCP = zeros(2,n);
j = 1;


for i=1:n
   
    % running simulation 
    [settings.CoeffsF,settings.CoeffsE,settings.Alphas,settings.Betas,settings.Altitudes,settings.Machs] = takefile(r_name{i});
    [apogee(i), max_a(i),Vexit(i),t,Xcp] = start_simulation(settings);
   
    % ottimizzazione 
    check=isnan(Xcp);
    while check(j)==1
        j=j+1;
    end
    
    XCP(i,1)=Xcp(j);
    XCP(i,2)=Xcp(end);
    
    
    % plot grafici stabilità 
    if grafico 
        plot(t,Xcp,'.')
        xlabel('time [s]')
        ylabel('S / M')
        hold on 
    end
end

% results control








