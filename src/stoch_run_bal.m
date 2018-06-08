function [LP,X,ApoTime] = stoch_run_bal(settings)
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

if not(settings.wind.model)
    if settings.wind.MagMin == settings.wind.MagMax && settings.wind.ElMin == settings.wind.ElMax
        error('In stochastic simulations the wind must setted with the random model, check config.m')
    end 
else
    if settings.wind.DayMin == settings.wind.DayMax && settings.wind.HourMin == settings.wind.HourMax
        error('In stochastic simulations with the wind model the day or the hour of launch must vary, check config.m')
    end
end

%% STARTING CONDITIONS

% Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

% PreAllocation
LP = zeros(settings.stoch.N,3);
X = zeros(settings.stoch.N,3);
ApoTime = zeros(settings.stoch.N,1);


%% PARFOR LOOP

parfor_progress(settings.stoch.N); % initiaize parfor loop
parpool;

parfor i = 1:settings.stoch.N
    
    %% WIND GENERATION
    
    if not(settings.wind.model)
        
    [uw,vw,ww] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
        settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
        settings.wind.MagMax);
    else
        Day = randi([settings.wind.DayMin,settings.wind.DayMax]);
        Hour = randi([settings.wind.HourMin,settings.wind.HourMax]);
        uw = 0; vw = 0; ww = 0;
    end

    %% ASCENT 

    [Ta,Ya] = ode113(@ascent,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
        settings,uw,vw,ww,Hour,Day);

    
    %% DESCEND
    
    % control if the second parachute fail
    if settings.sdf
        
        % first parachute descent phase
        
        para = 1; % Flag for Drogue 1
        X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
        [~,Yd1] = ode113(@descent_parachute,settings.ode.timedrg1,X0d1,...
            settings.ode.optionsdrg1,settings,uw,vw,ww,para);
        
        % after failure of drogue 2 ballistic descent
        
        Q0 = angle2quat(90*pi/180,0,0,'ZYX')';
        X0b = [Yd1(end,:) 0 0 0 Q0'];
        [~,Yb] = ode113(@descent_ballistic,settings.ode.timedesc,X0b,settings.ode.optionsdesc,...
            settings,uw,vw,ww,Hour,Day);
        
    else
        % total ballistic descend, so no drogue will be used
        
        [~,Yd] = ode113(@descent_ballistic,settings.ode.timedesc,Ya(end,1:13),settings.ode.optionsdesc,...
            settings,uw,vw,ww,Hour,Day);
    end


    %% FINAL STATE ASSEMBLING
    
    %Total State
    if settings.sdf
        LP(i,:) = Yb(end,1:3);
    else
        LP(i,:) = Yd(end,1:3);
    end
    
    
    X(i,:) = [Ya(end,1); Ya(end,2); -Ya(end,3)]
    ApoTime(i) = Ta(end);
    
    parfor_progress;

end

