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

%% STARTING CONDITIONS

% Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

%% WIND GENERATION

    [uw,vw,ww,~] = wind_const_generator(settings.wind.AzMin,settings.wind.AzMax,...
    settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
    settings.wind.MagMax);

tf = settings.ode.final_time;

%% ASCENT

[Ta,Ya] = ode113(@ascent,[0,tf],X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww);
[data_ascent] = RecallOdeFcn(@ascent,Ta,Ya,settings,uw,vw,ww);
data_ascent.state.Y = Ya;
data_ascent.state.T = Ta;
save('data_ascent.mat', 'data_ascent');

%% FINAL STATE ASSEMBLING 
    
% Total State
Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6)) Ya(:,7:13)];
% Total Time
Tf = [Ta];

%% TIME, POSITION AND VELOCITY AT APOGEE

bound_value.td1 = Ta(end);
bound_value.Xd1 = [Ya(end,2), Ya(end,1), -Ya(end,3)];
bound_value.Vd1 = quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6));

