function [p,flag] = LaunchProb(settings,data_ascent,data_descent)
% LaunchProb - This function allows to compute the probability of launch
% subject to the constraints specified in config.m
%
% CONSTRAINTS:  
%               - Stability Margin                                                     #1    
%               - Landing area                                                         #2     
%               - Crash of 1st drogue at opening (apogee + delay time of opening)      #3     
%     
% OUTPUTS
%               - p: probability
%               - flag: outputs if each constrain is active as follows
%
% FLAGS (logical explanation):     #1  #2  #3
%                                  1   1   1    ==> flag = 1
%                                  1   1   0    ==> flag = -1
%                                  1   0   1    ==> flag = -2
%                                  0   1   1    ==> flag = -3
%                                  1   0   0    ==> flag = -4
%                                  0   1   0    ==> flag = -5
%                                  0   0   1    ==> flag = -6
%                                  0   0   0    ==> flag = -7

% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
% Release date: 16/11/2018

N_simul = settings.stoch.N;
x_lim = - settings.stoch.x_lim;
% t_delay = settings.stoch.prob.t_delay;
P_lim = settings.stoch.P_lim;

flag = zeros(N_simul,1); 

for i = 1:N_simul
%     N_asc = length(data_ascent{i}.integration.t);
    
    flag_XCP = all(-data_ascent{i}.coeff.XCP(not(isnan(data_ascent{i}.coeff.XCP))) > 0);
    flag_land = data_descent{i}.state(3).Y(end,2) > x_lim;
    
    V = norm(data_ascent{i}.state.Y(end,4:6) - data_ascent{i}.wind.body_wind(1:3,end));
    P_apo = 0.5*data_ascent{i}.air.rho(end)*V^2*settings.para1.S*settings.para1.CD;
    flag_d1 = (P_apo < P_lim);
    
    if flag_XCP == 1
        
        if flag_land == 1
            
            if flag_d1 == 1
                flag(i) = 1;
            else 
                flag(i) = -1;
            end
            
        else
            
            if flag_d1 == 1
                flag(i) = -2;
            else
                flag(i) = -4;
            end
            
        end
       
    else
        
        if flag_land == 1
            
            if flag_d1 == 1
                flag(i) = -3;
            else
                flag(i) = -5;
            end
            
        else
            
            if flag_d1 == 1
                flag(i) = -6;
            else
                flag(i) = -7;
            end
            
        end
        
    end
 
end
    
k = length(flag(flag == 1));
p = k/N_simul;
