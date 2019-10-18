function [apogee,t,vect_XCP] = run_sim(settings)

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


%% CALCULATE OUTPUT QUANTITIES 

[apogee,] = max(-Ya(:,3));
t=Ta;
vect_XCP=-data_ascent.coeff.XCP;
 

