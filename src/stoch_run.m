function [LP,X,ApoTime,data_ascent,data_para] = stoch_run(settings)
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

if settings.OMEGAmin == settings.OMEGAmax && settings.PHImin == settings.PHImax
    if not(settings.wind.model) && not(settings.wind.input)
        
        if settings.wind.MagMin == settings.wind.MagMax && settings.wind.ElMin == settings.wind.ElMax
            error('In stochastic simulations the wind must setted with the random model, check config.m')
        end
        
    elseif settings.wind.model
        
        if settings.wind.DayMin == settings.wind.DayMax && settings.wind.HourMin == settings.wind.HourMax
            error('In stochastic simulations with the wind model the day or the hour of launch must vary, check config.m')
        end
        
    end
    
    if settings.wind.input && all(settings.wind.input_uncertainty) == 0
        error('In stochastic simulations the wind input model, the uncertainty must be different to 0 check config.m')
    end
end

%% STARTING CONDITIONS

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';

%PreAllocation
LP = zeros(settings.stoch.N,3);
X = zeros(settings.stoch.N,3);
ApoTime = zeros(settings.stoch.N,1);

tf = settings.ode.final_time;

%% PARFOR LOOP
parfor_progress(settings.stoch.N);
parpool;

parfor i = 1:settings.stoch.N
    
    %% WIND GENERATION
    
    if not(settings.wind.model)
        
        Day = 0; Hour = 0;
        
        if settings.wind.input
            if settings.wind.input_uncertainty == 0
                uncert = [0,0];
            else
                
                signn = randi([0,1]);
                
                if signn
                    unc = - settings.wind.input_uncertainty;
                else
                    unc = settings.wind.input_uncertainty;
                end
                
                 uncert = rand(1,2).*unc;
            end
            uw = 0; vw = 0; ww = 0;
        else
            [uw,vw,ww,Azw] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
                settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
                settings.wind.MagMax);
            uncert = [0; 0];
        end
        
    else
        Day = randi([settings.wind.DayMin,settings.wind.DayMax]);
        Hour = randi([settings.wind.HourMin,settings.wind.HourMax]);
        uw = 0; vw = 0; ww = 0; uncert = [0,0];
    end
    
    
    %% ASCENT
    % Attitude
    OMEGA = settings.OMEGAmin + rand*(settings.OMEGAmax - settings.OMEGAmin);
    PHI = settings.PHImin + rand*(settings.PHImax - settings.PHImin);
    if settings.upwind
        PHI = Azw + 180;
    end
    
    Q0 = angle2quat(PHI,OMEGA,0*pi/180,'ZYX')';
    X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];
    
    [Ta,Ya] = ode113(@ascent,[0,tf],X0a,settings.ode.optionsasc,...
        settings,uw,vw,ww,uncert,Hour,Day,OMEGA);
    [data_ascent{i}] = RecallOdeFcn(@ascent,Ta,Ya,settings,uw,vw,ww,uncert,Hour,Day,OMEGA);
    data_ascent{i}.state.Y = Ya;
    data_ascent{i}.state.T = Ta;
    
    if not(settings.ao)
        
        %% DROGUE 1
        % Initial Condition are the last from ascent (need to rotate because
        % velocities are outputted in body axes)
        
        para = 1; % Flag for Drogue 1
        X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
        
        if settings.rocket_name == "R2A_hermes" && settings.ldf
            [Td1,Yd1] = ode113(@descent_parachute,[Ta(end),tf],X0d1,settings.ode.optionsdrg2,settings,uw,vw,ww,para,uncert);
        else
            [Td1,Yd1] = ode113(@descent_parachute,[Ta(end),tf],X0d1,settings.ode.optionsdrg1,settings,uw,vw,ww,para,uncert);
        end
        
        [data_para1] = RecallOdeFcn(@descent_parachute,Td1,Yd1,settings,uw,vw,ww,para,uncert);
        data_para1.state.Y = Yd1;
        data_para1.state.T = Td1;
        
        %% DROGUE 2
        
        para = 2; % Flag for Drogue 2
        X0d2 = Yd1(end,:); % Initial Condition are the last from drogue 1 descent
        
        if settings.rocket_name == "R2A" || (settings.rocket_name == "R2A_hermes"  && not(settings.ldf))
            [Td2,Yd2] = ode113(@descent_parachute,[Td1(end),tf],X0d2,...
                settings.ode.optionsdrg2,settings,uw,vw,ww,para,uncert);
            [data_para2] = RecallOdeFcn(@descent_parachute,Td2,Yd2,settings,uw,vw,ww,para,uncert);
            data_para2.state.Y = Yd2;
            data_para2.state.T = Td2;
        end
        
        %% ROGALLO WING
        % Initial Condition are the last from drogue 2 descent
        
        if settings.rocket_name == "R2A"
            
            if not(settings.ldf)
                para = 3;             % Flag for Main (Rogall)
            end
            
            X0m = Yd2(end,:);
            [T_rog,Y_rog] = ode113(@descent_parachute,[Td2(end),tf],X0m,...
                settings.ode.optionsrog,settings,uw,vw,ww,para,uncert,Hour,Day);
            [data_para3] = RecallOdeFcn(@descent_parachute,T_rog,Y_rog,settings,uw,vw,ww,para,uncert,Hour,Day);
            data_para3.state.Y = Y_rog;
            data_para3.state.T = T_rog;
            data_para{i} = cell2struct(cellfun(@vertcat,struct2cell(data_para1),struct2cell(data_para2),struct2cell(data_para3),'uni',0),fieldnames(data_para1),1);
        end
    end
    
    %% FINAL STATE ASSEMBLING
    
    switch settings.rocket_name
        
        case 'R2A'
            if not(settings.ao)
                Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2;Y_rog];
                LP(i,:) = Yf(end,1:3);
            end
            
        case 'R2A_hermes'
            
            % Total State
            if not(settings.ao)
                if not(settings.ldf)
                    Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2];
                    data_para{i} = cell2struct(cellfun(@vertcat,struct2cell(data_para1),struct2cell(data_para2),'uni',0),fieldnames(data_para1),1);
                else
                    Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1];
                     data_para{i} = data_para1;
                end
                LP(i,:) = Yf(end,1:3);
            end
            
            X(i,:) = [Ya(end,1); Ya(end,2); -Ya(end,3)]
            ApoTime(i) = Ta(end);
            
            parfor_progress;
            
    end
    
end

end
