clc; clear; close all
run('config_R2A_hermes.m')

heights = 200:100:1200;
magn_w = 1:10;

Nh = length(heights);
Nw = length(magn_w);
R = zeros(Nh,Nw);

for i = 1:Nh
    
    for j = 1:Nw
        settings.wind.MagMin = magn_w(j);
        settings.wind.MagMax = magn_w(j);
        settings.zdrg2 = heights(i);
        [Tf,Yf,Ta,Ya,bound_value] = std_run(settings);
        R(i,j) = norm([Yf(end,1),Yf(end,2)]);
        
    end
end

R_round = round(R,2,'significant');
varNames = strcat('h_Rog_', string(200:100:1200));
varNames = arrayfun(@(x)char(varNames(x)),1:numel(varNames),'uni',false);
RowNames = strcat('w_ ', string(1:10));
RowNames = arrayfun(@(x)char(RowNames(x)),1:numel(RowNames),'uni',false);
T = array2table(R_round','RowNames',RowNames,'VariableNames',varNames);