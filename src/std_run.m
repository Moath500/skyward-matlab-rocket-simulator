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
if settings.sdf
    warning('The second drogue failure can be simulated just in ballistic simulations, check settings.sdf & settings.ballistic in config.m')
end

if settings.wind.HourMin ~= settings.wind.HourMax || settings.wind.HourMin ~= settings.wind.HourMax
    error('In standard simulations with the wind model the day and the hour of launch must be unique, check config.m')
end

%% STARTING CONDITIONS

% Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

%% WIND GENERATION

if settings.wind.model || settings.wind.input   % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0;
else
    [uw,vw,ww] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
        settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
        settings.wind.MagMax);
    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
end

%% ASCENT

[Ta,Ya] = ode113(@ascent,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww);

%% DROGUE 1
% Initial Condition are the last from ascent (need to rotate because
% velocities are outputted in body axes)

para = 1; % Flag for Drogue 1
X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
[Td1,Yd1] = ode113(@descent_parachute,settings.ode.timedrg1,X0d1,...
    settings.ode.optionsdrg1,settings,uw,vw,ww,para);

%% DROGUE 2
% Initial Condition are the last from drogue 1 descent

para = 2; % Flag for Drogue 2
X0d2 = Yd1(end,:);
[Td2,Yd2] = ode113(@descent_parachute,settings.ode.timedrg2,X0d2,...
    settings.ode.optionsdrg2,settings,uw,vw,ww,para);

%% ROGALLO WING
% Initial Condition are the last from drogue 2 descent

if not(settings.ldf)
para = 3;             % Flag for Main (Rogall)
end

X0m = Yd2(end,:);
[Trog,Yrog] = ode113(@descent_parachute,settings.ode.timerog,X0m,...
    settings.ode.optionsrog,settings,uw,vw,ww,para);



%% FINAL STATE ASSEMBLING

% Total State
Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yd2;Yrog];

% Total Time
Tf = [Ta; Ta(end)+Td1; Ta(end)+Td1(end)+Td2;Ta(end)+Td1(end)+Td2(end)+Trog];

%% TIME, POSITION AND VELOCITY AT DROGUES DEPLOYMENT

bound_value.td1 = Ta(end);
bound_value.td2 = Ta(end)+Td1(end);
bound_value.trog = Ta(end)+Td1(end)+Td2(end);
bound_value.Xd1 = [Ya(end,2), Ya(end,1), -Ya(end,3)];
bound_value.Xd2 = [Yd1(end,2), Yd1(end,1), -Yd1(end,3)];
bound_value.Xrog = [Yd2(end,2), Yd2(end,1), -Yd2(end,3)];
bound_value.Vd1 = quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6));
bound_value.Vd2 = [Yd1(end,4) Yd1(end,5) -Yd1(end,6)];
bound_value.Vrog = [Yd2(end,4) Yd2(end,5) -Yd2(end,6)];

