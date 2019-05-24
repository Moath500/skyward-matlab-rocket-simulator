function [p,flag,ind_Pin,ind_Pout,RP] = LaunchProb(settings,data_ascent,data_descent,LP)
%{
LaunchProb - This function allows to compute the probability of launch
subject to the constraints specified in config.m

CONSTRAINTS:
              - Stability Margin                                                     #1
              - Landing area                                                         #2
              - Crash of 1st drogue at opening (apogee + delay time of opening)      #3

OUTPUTS
              - p: probability
              - flag: outputs if each constrain is active as follows

FLAGS (logical explanation):     #1  #2  #3
                                 1   1   1    ==> flag = 1
                                 1   1   0    ==> flag = -1
                                 1   0   1    ==> flag = -2
                                 0   1   1    ==> flag = -3
                                 1   0   0    ==> flag = -4
                                 0   1   0    ==> flag = -5
                                 0   0   1    ==> flag = -6
                                 0   0   0    ==> flag = -7

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 16/11/2018

%}

N_simul = settings.stoch.N;
V_lim = settings.stoch.prob.V_lim;

flag = zeros(N_simul,1);

if settings.rocket_name == "R2A_hermes"
    % pre-allocation
    dist = zeros(settings.stoch.N,1);
    RP = zeros(settings.stoch.N,2);
    
    a = settings.prob.SafeEllipse.a;
    b = settings.prob.SafeEllipse.b;
    x0 = settings.prob.SafeEllipse.x0;
    y0 = settings.prob.SafeEllipse.y0;
    alpha = settings.prob.SafeEllipse.alpha;
    c = sqrt(b^2 - a^2);                        % [m] safety ellipse focii distance
    F1 = [x0, y0+c];                            % [m] first safety ellipse not-rotated foci coordinates
    F2 = [x0, y0-c];                            % [m] second safety ellipse not-rotated foci coordinates
    
    R = [cosd(alpha), - sind(alpha)             % Rotation matrix
        sind(alpha),   cosd(alpha)];
    
    F1_hat = R*F1';                             % [m] first rotated focii
    F2_hat = R*F2';                             % [m] second rotated focii
    
else
    x_lim = - settings.stoch.prob.x_lim;
end

for i = 1:N_simul
    
    flag_XCP = all(-data_ascent{i}.coeff.XCP(not(isnan(data_ascent{i}.coeff.XCP))) > 0);
    
    if settings.rocket_name == "R2A"
        flag_land = data_descent{i}.state(3).Y(end,2) > x_lim;
    else
        
        P = LP(:,[2,1]);
        dist(i) = norm(P(i,:) - F1_hat') + norm(P(i,:) - F2_hat');
        flag_land = dist(i) < 2*b;
        ind_Pin = find(dist < 2*b);
        ind_Pout = find(dist > 2*b);
        RP(i,:) = data_descent{1, i}.state(1).Y(end,1:2);
        
    end
    
    %     t_apogee = data_ascent{i}.integration.t(end);
    %     [~, ind_ComputedApogee] = min(abs(data_descent{1, i}.integration(1).t(:) - t_apogee - t_delay));
    %     V_apo = norm(data_descent{i}.state(1).Y(ind_ComputedApogee,4:6) - ...
    %         data_descent{i}.wind(1).body_wind(1:3,ind_ComputedApogee));
    
    V_apo = norm(data_ascent{i}.state.Y(end,4:6) - ...
            data_ascent{i}.wind.body_wind(1:3,end)');
    
    flag_drogue1 = V_apo < V_lim;
    
    if flag_XCP == 1
        
        if flag_land == 1
            
            if flag_drogue1 == 1
                flag(i) = 1;
            else
                flag(i) = -1;
            end
            
        else
            
            if flag_drogue1 == 1
                flag(i) = -2;
            else
                flag(i) = -4;
            end
            
        end
        
    else
        
        if flag_land == 1
            
            if flag_drogue1 == 1
                flag(i) = -3;
            else
                flag(i) = -5;
            end
            
        else
            
            if flag_drogue1 == 1
                flag(i) = -6;
            else
                flag(i) = -7;
            end
            
        end
        
    end
    
end

k = length(flag(flag == 1));
p = k/N_simul*100;
