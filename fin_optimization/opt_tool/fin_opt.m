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

run Auto_Matrices.m


run config.m

n=length(data);

% settings 
grafico = true;
cal_min = 1;                    % minum stability margin required

%% RUN 

% preallocation 
apogee = zeros(1,n);
XCP = zeros(2,n);
leg = cell(1,n);
j = 1;

tic 
for i=1:n
   
    % running simulation 
    [settings.CoeffsF,settings.CoeffsE,settings.Alphas,settings.Betas,settings.Altitudes,settings.Machs] = takefile(data{i});
    [apogee(i),t,Xcp] = run_sim(settings);
   
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
        
        leg{i}=strcat('Matrice_',int2str(i));
    end
end
legend(leg);

%% results control

index=find(XCP(1:n,1)>cal_min);
v_xcp=XCP(index,1);

apo_ok=apogee(index);
[apo_max,ind_max]=max(apo_ok);

ind = index(ind_max); 

%% print results 

fprintf('il massimo apogeo è %g ottenuto con le alette: \n',apo_max)
fprintf('forma a %s \n',data{ind}.shape)
fprintf('corda max: %g [m] \n',data{ind}.c_max)
fprintf('corda min: %g [m] \n',data{ind}.c_min)
fprintf('altezza: %g [m] \n',data{ind}.h)


toc 




