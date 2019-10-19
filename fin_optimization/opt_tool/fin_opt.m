%{

FIN_OPT - This script chose the best fin set between the chosen.
For every set, missile datcom is used to predict the aerodynamic coefficients.

Author: Matteo Pozzoli
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: matteo.pozzoli@skywarder.eu
Website: http://www.skywarder.eu
Release date: 14/10/2019

%}

clear 
close all
clc 

path = genpath(pwd);
addpath(path);


%% COMPUTING AERODYNAMIC DATA
run Auto_Matrices.m

%% RETRIVING GEOMETRICAL DATA
run config.m

%% FINSETs SIMULATION 
tic

% preallocation 
n = length(data);
apogee = zeros(1, n);
XCP = zeros(2, n);

figure; hold on; xlabel('time [s]'); ylabel('S / M');

for i = 1:n
    data_ith = data{i};
    [settings] = ManageData(data_ith, settings);
    [apogee(i), t, Xcp] = run_sim(settings);
   
    % optimization values 
    ind_notNan = (not(isnan(Xcp)));
    Xcp = Xcp(ind_notNan);
    t = t(ind_notNan);
    XCP(i, 1) = Xcp(1);
    XCP(i, 2) = Xcp(end);
    
    % plot grafici stabilità 
    if settings.plot
        plot(t, Xcp, '.')
    end
end

legend(strcat('Matrix_', string(1:n)));

%% OPTIMIZATION
% chosing between XPC and apogee
indexes = find(XCP(1:n, 1) > settings.cal_min);
v_xcp = XCP(indexes, 1);
apo_ok = apogee(indexes);
[apo_max, ind_max] = max(apo_ok);

best = indexes(ind_max); 

FOtime = toc;

%% print results 
fprintf('COMPUTATIONAL EFFORT: \n\n')
fprintf('- AutoMatrices time, %g [s]\n', AMtime)
fprintf('- Simulations time, %g [s]\n', FOtime)
fprintf('- Total time, %g [s]\n\n\n', FOtime + AMtime)
fprintf('BEST FINSET: \n\n')
fprintf('- shape, %s \n', data{best}.shape)
fprintf('- attached chord, %g [m] \n', data{best}.c_max)
fprintf('- free chord, %g [m] \n', data{best}.c_min)
fprintf('- height, %g [m] \n\n\n', data{best}.h)
fprintf('BEST FINSET SIMULATION RESULTS: \n\n')
fprintf('- apogee, %g [m]: \n', apo_max)
fprintf('- stability margin @launchpad exit, %g \n', XCP(best, 1))
fprintf('- stability margin @apogee, %g \n', XCP(best, 2))
