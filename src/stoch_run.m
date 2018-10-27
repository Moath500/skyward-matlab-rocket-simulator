function [LP,X,ApoTime] = stoch_run(settings)
%STD RUN - This function runs a stochastic simulation (parallel)
% OUTPUTS
% LP: Landing Points
% Z: Apogee Altitudes

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 29.V.2014
% License:  2-clause BSD

warning off

if settings.wind.model && settings.wind.input
    error('Both wind model and input wind are true, select just one of them')
end

if not(settings.wind.model) && not(settings.wind.input)
    
    if settings.wind.MagMin == settings.wind.MagMax && settings.wind.ElMin == settings.wind.ElMax
        error('In stochastic simulations the wind must setted with the random model, check config.m')
    end
    
elseif settings.wind.model
    
    if settings.wind.DayMin == settings.wind.DayMax && settings.wind.HourMin == settings.wind.HourMax
        error('In stochastic simulations with the wind model the day or the hour of launch must vary, check config.m')
    end
    
end

if settings.wind.input && settings.wind.input_uncertainty == 0
    error('In stochastic simulations the wind input model, the uncertainty must be different to 0 check config.m')
end

%% STARTING CONDITIONS

% Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

%PreAllocation
LP = zeros(settings.stoch.N,3);
X = zeros(settings.stoch.N,3);
ApoTime = zeros(settings.stoch.N,1);

%% PARFOR LOOP
parfor_progress(settings.stoch.N);
parpool;

parfor i = 1:settings.stoch.N
    
    %% WIND GENERATION
    
    if not(settings.wind.model)
        
        Day = 0; Hour = 0;
        
        if settings.wind.input
            uncert = randi(settings.wind.input_uncertainty,[1,3]);
            uw = 0; vw = 0; ww = 0;
        else
            [uw,vw,ww] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
                settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
                settings.wind.MagMax);
            uncert = 0;
        end
        
    else
        Day = randi([settings.wind.DayMin,settings.wind.DayMax]);
        Hour = randi([settings.wind.HourMin,settings.wind.HourMax]);
        uw = 0; vw = 0; ww = 0; uncert = [0,0];
    end
    

    
    %% ASCENT

    [Ta,Ya] = ode113(@ascent,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
        settings,uw,vw,ww,uncert,Hour,Day);

    %% DROGUE 1
    % Initial Condition are the last from ascent (need to rotate because
    % velocities are outputted in body axes)

    para = 1; % Flag for Drogue 1
    X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
    [~,Yd1] = ode113(@descent_parachute,settings.ode.timedrg1,X0d1,...
        settings.ode.optionsdrg1,settings,uw,vw,ww,para,uncert,Hour,Day);

    %% DROGUE 2 
    % Initial Condition are the last from drogue 1 descent
    
    para = 2; %Flag for Drogue 2
    X0d2 = Yd1(end,:);
    [~,Yd2] = ode113(@descent_parachute,settings.ode.timedrg2,X0d2,...
        settings.ode.optionsdrg2,settings,uw,vw,ww,para,uncert,Hour,Day);

    %% ROGALLO WING
    % Initial Condition are the last from drogue 2 descent
    
    switch settings.rocket_name
        
        case 'R2A'
            
            if not(settings.ldf)
                para = 3;             % Flag for Main (Rogall)
            end
            
            X0m = Yd2(end,:);
            [~,Ym] = ode113(@descent_parachute,settings.ode.timerog,X0m,...
                settings.ode.optionsrog,settings,uw,vw,ww,para,uncert,Hour,Day);
            
            %% FINAL STATE ASSEMBLING
            
            if not(settings.ao)
                Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1; Yd2;Ym];
                LP(i,:) = Yf(end,1:3);
            end
            
        case 'R2A_hermes'
            
            % Total State
            if not(settings.ao)
                Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2];
                LP(i,:) = Yf(end,1:3);
            end
            
    end
    
    X(i,:) = [Ya(end,1); Ya(end,2); -Ya(end,3)]
    ApoTime(i) = Ta(end);

    parfor_progress;

end
end
