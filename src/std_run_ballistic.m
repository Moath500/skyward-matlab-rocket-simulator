function [Tf,Yf,Ta,Ya,bound_value] = std_run_ballistic(settings)
% STD RUN BALLISTIC - This function runs a ballistic (non-stochastic) simulation
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

if settings.rocket_name == 'R2A'
    if settings.ldf
        error('Landing with the second drogue can be simulated just in standard simulations, check settings.ldf & settings.ballistic in config.m')
    end
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

if settings.wind.input && settings.wind.input_uncertainty ~= 0
    uncert = randi(settings.wind.input_uncertainty,[1,2]);
else
    uncert = [0,0];
end
    

%% ASCENT

[Ta,Ya] = ode113(@ascent,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww,uncert);

%% DESCEND 

% control if the second parachute fail
if settings.sdf
    % first parachute descent phase 
    
    para = 1; % Flag for Drogue 1
    X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
    [Td1,Yd1] = ode113(@descent_parachute,settings.ode.timedrg1,X0d1,...
    settings.ode.optionsdrg1,settings,uw,vw,ww,para,uncert);

    % ballistic descent after failure of drogue 2 
    
    X0b1 = Yd1(end,:);
    Q0 = angle2quat(0,0,0,'ZYX')';
    X0b = [X0b1,0,0,0,Q0'];
    [Tb,Yb] = ode45(@descent_ballistic,settings.ode.timedesc,X0b,settings.ode.optionsdesc,...
        settings,uw,vw,ww);
    
else 

% total ballistic descend, so no drogue will be used

[Td,Yd] = ode113(@descent_ballistic,settings.ode.timedesc,Ya(end,1:13),settings.ode.optionsdesc,...
    settings,uw,vw,ww,uncert);
end

%% FINAL STATE ASSEMBLING 

if settings.sdf
% Total State
Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1;Yb(:,1:3),...
      quatrotate(quatconj(Yb(:,10:13)),Yb(:,4:6))];
% Total Time
Tf = [Ta; Ta(end)+Td1; Ta(end)+Td1(end)+Tb];

else
    
% Total State
Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6)) Ya(:,7:13)
      Yd(:,1:3) quatrotate(quatconj(Yd(:,10:13)),Yd(:,4:6)) Yd(:,7:13)];
% Total Time
Tf = [Ta; Ta(end)+Td];

end

%% TIME, POSITION AND VELOCITY AT APOGEE

bound_value.td1 = Ta(end);
bound_value.Xd1 = [Ya(end,2), Ya(end,1), -Ya(end,3)];
bound_value.Vd1 = quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6));

if settings.sdf 
    bound_value.td2 = Ta(end)+Td1(end);
    bound_value.Xd2 = [Yd1(end,2), Yd1(end,1), -Yd1(end,3)];
    bound_value.Vd2 = [Yd1(end,4), Yd1(end,5), -Yd1(end,6)];
end

