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

[Ta,Ya] = ode113(@ascent,[0,tf],X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww,uncert);
[data_ascent] = RecallOdeFcn(@ascent,Ta,Ya,settings,uw,vw,ww,uncert);
data_ascent.state.Y = Ya;
data_ascent.state.T = Ta;
save('ascent_plot.mat', 'data_ascent');

%% DESCEND 

[Td,Yd] = ode113(@descent_ballistic,[Ta(end),tf],Ya(end,1:13),settings.ode.optionsdesc,...
    settings,uw,vw,ww,uncert);
[data_bal] = RecallOdeFcn(@descent_ballistic,Td,Yd,settings,uw,vw,ww,uncert);
data_bal.state.Y = Yd;
data_bal.state.T = Td;
save('descent_plot.mat', 'data_bal');

%% FINAL STATE ASSEMBLING 
    
% Total State
Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6)) Ya(:,7:13)
      Yd(:,1:3) quatrotate(quatconj(Yd(:,10:13)),Yd(:,4:6)) Yd(:,7:13)];
% Total Time
Tf = [Ta; Ta(end)+Td];

%% TIME, POSITION AND VELOCITY AT APOGEE

bound_value.td1 = Ta(end);
bound_value.Xd1 = [Ya(end,2), Ya(end,1), -Ya(end,3)];
bound_value.Vd1 = quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6));

