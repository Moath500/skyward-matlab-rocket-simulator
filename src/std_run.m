function [Tf,Yf,Ta,Ya,bound_value] = std_run(settings)
% STD RUN - This function runs a standard (non-stochastic) simulation
% OUTPUTS
% Tf: Time steps
% Yf: Final State

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 31.XII.2014
% License:  2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

if settings.wind.model && settings.wind.input
    error('Both wind model and input wind are true, select just one of them')
end

if settings.wind.HourMin ~= settings.wind.HourMax || settings.wind.HourMin ~= settings.wind.HourMax
    error('In standard simulations with the wind model the day and the hour of launch must be unique, check config.m')
end

if settings.OMEGAmin ~= settings.OMEGAmax || settings.PHImin ~= settings.PHImax 
    error('In a single simulation the launchpad configuration has to be unique, check config.m')
end

%% WIND GENERATION
if settings.wind.model || settings.wind.input   % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0; 
else
    [uw,vw,ww,Azw] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
        settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
        settings.wind.MagMax);
    
    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
    
end

if settings.wind.input && all(settings.wind.input_uncertainty) ~= 0
    signn = randi([1,4]); % 4 sign cases
    unc = settings.wind.input_uncertainty;
    
    switch signn
        case 1
            %                       unc = unc;
        case 2
            unc(1) = - unc(1);
        case 3
            unc(2) = - unc(2);
        case 4
            unc = - unc;
    end
    
    uncert = rand(1,2).*unc;
else
    uncert = [0,0];
end

tf = settings.ode.final_time;


%% STARTING CONDITIONS
% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';

settings.OMEGA = settings.OMEGAmin;

Azw
% Attitude
if settings.wind.input || settings.wind.model
    settings.PHI = settings.PHImin;
else
    
    if settings.upwind
        settings.PHI = mod(Azw + pi, 2*pi);
    else
        settings.PHI = settings.PHImin + rand*(settings.PHImax - settings.PHImin);
    end
    
end

Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

%% ASCENT
% checking if the actuation delay is different from zero
if settings.para1.delay ~= 0
    [Ta1,Ya1] = ode113(@ascent,[0,tf],X0a,settings.ode.optionsasc1,settings,uw,vw,ww,uncert);
    
    [Ta2,Ya2] = ode113(@ascent,[Ta1(end),Ta1(end) + settings.para1.delay],Ya1(end,:),settings.ode.optionsasc2,settings,uw,vw,ww,uncert);
    Ta = [Ta1; Ta2(2:end)];
    Ya = [Ya1; Ya2(2:end,:)];
else
    [Ta,Ya] = ode113(@ascent,[0,tf],X0a,settings.ode.optionsasc,settings,uw,vw,ww,uncert);
end

[data_ascent] = RecallOdeFcn(@ascent,Ta,Ya,settings,uw,vw,ww,uncert);
data_ascent.state.Y = Ya;
data_ascent.state.T = Ta;
save('ascent_plot.mat', 'data_ascent');

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
% Initial Condition are the last from drogue 1 descent

para = 2; % Flag for Drogue 2
X0d2 = Yd1(end,:);

if settings.rocket_name == "R2A" || (settings.rocket_name == "R2A_hermes" && not(settings.ldf))
    [Td2,Yd2] = ode113(@descent_parachute,[Td1(end),tf],X0d2,...
        settings.ode.optionsdrg2,settings,uw,vw,ww,para,uncert);
    [data_para2] = RecallOdeFcn(@descent_parachute,Td2,Yd2,settings,uw,vw,ww,para,uncert);
    data_para2.state.Y = Yd2;
    data_para2.state.T = Td2;
end

%% ROGALLO WING AND FINAL STATE ASSEMBLING
% Initial Condition are the last from drogue 2 descent

switch settings.rocket_name
    
    case 'R2A'
        
        if not(settings.ldf)
            para = 3;             % Flag for Main (Rogall)
        end
        
        X0m = Yd2(end,:);
        [Trog,Yrog] = ode113(@descent_parachute,[Td2(end),tf],X0m,settings.ode.optionsrog,settings,uw,vw,ww,para,uncert);
        [data_para3] = RecallOdeFcn(@descent_parachute,Trog,Yrog,settings,uw,vw,ww,para,uncert);
        data_para3.state.Y = Yrog;
        data_para3.state.T = Trog;

        % Total State
        Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2;Yrog];
        
        % Total Time
        Tf = [Ta; Td1; Td2; Trog];
        
        data_para = cell2struct(cellfun(@vertcat,struct2cell(data_para1),struct2cell(data_para2),struct2cell(data_para3),'uni',0),fieldnames(data_para1),1);
        
    case 'R2A_hermes'
        
        if not(settings.ldf)
            % Total State
            Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2];
            % Total Time
            Tf = [Ta; Td1; Td2];
            data_para = cell2struct(cellfun(@vertcat,struct2cell(data_para1),struct2cell(data_para2),'uni',0),fieldnames(data_para1),1);
        else 
            % Total State
            Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1];
            % Total Time
            Tf = [Ta; Td1];
            data_para = data_para1;
        end
end

save('descent_para_plot.mat', 'data_para')

%% TIME, POSITION AND VELOCITY AT DROGUES DEPLOYMENT

bound_value.td1 = Ta(end);
bound_value.td2 = Td1(end);
bound_value.Xd1 = [Ya(end,2), Ya(end,1), -Ya(end,3)];
bound_value.Xd2 = [Yd1(end,2), Yd1(end,1), -Yd1(end,3)];
bound_value.Vd1 = quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6));
bound_value.Vd2 = [Yd1(end,4) Yd1(end,5) -Yd1(end,6)];

if settings.rocket_name == "R2A_hermes" && not(settings.ldf) || settings.rocket_name == "R2A" && not(settings.ldf)
    bound_value.Xrog = [Yd2(end,2), Yd2(end,1), -Yd2(end,3)];
    bound_value.trog = Td2(end);
    bound_value.Vrog = [Yd2(end,4) Yd2(end,5) -Yd2(end,6)];
end



