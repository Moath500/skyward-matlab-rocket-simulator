function [Tf, Yf, Ta, Ya, bound_value] = std_run_ballistic(settings)
%{

 STD RUN BALLISTIC - This function runs a ballistic (non-stochastic) simulation

OUTPUTS:
            - Ta: Ascent time steps
            - Ya: Ascent State

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept
email: adriano.filippo.inno@skywarder.eu
Revision date: 09/10/2019

%}

if settings.wind.model && settings.wind.input
    error('Both wind model and input wind are true, select just one of them')
end

if settings.wind.HourMin ~= settings.wind.HourMax || settings.wind.HourMin ~= settings.wind.HourMax
    error('In standard simulations with the wind model the day and the hour of launch must be unique, check config.m')
end

if settings.OMEGAmin ~= settings.OMEGAmax || settings.PHImin ~= settings.PHImax 
    error('In a single simulation the launchpad configuration has to be unique, check config.m')
end

%% STARTING CONDITIONS

% Attitude
settings.OMEGA = settings.OMEGAmin;
settings.PHI = settings.PHImin;

if settings.upwind
    settings.PHI = mod(Azw + pi, 2*pi);
end

% Attitude
Q0 = angle2quat(settings.PHI, settings.OMEGA, 0*pi/180, 'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
Y0a = [X0; V0; W0; Q0; settings.m0; settings.Ixxf; settings.Iyyf; settings.Izzf];

%% WIND GENERATION

if settings.wind.model || settings.wind.input   % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0;
else 
    [uw,vw,ww,~] = wind_const_generator(settings.wind.AzMin, settings.wind.AzMax,...
    settings.wind.ElMin, settings.wind.ElMax, settings.wind.MagMin, settings.wind.MagMax);

    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
    
end

if settings.wind.input && all(settings.wind.input_uncertainty ~= 0)
    signn = randi([1, 4]); % 4 sign cases
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

%% ASCENT
% ascent phase computation
[Ta, Ya] = ode113(@ascent, [0, tf], Y0a, settings.ode.optionsasc1, settings, uw, vw, ww, uncert);
[data_ascent] = RecallOdeFcn(@ascent, Ta, Ya, settings, uw, vw, ww, uncert);
data_ascent.state.Y = Ya;
data_ascent.state.T = Ta;
save('ascent_plot.mat', 'data_ascent');

%% DESCEND 
% Initial Condition are the last from ascent
[Td,Yd] = ode113(@descent_ballistic, [Ta(end), tf], Ya(end, 1:13), settings.ode.optionsdesc, settings, uw, vw, ww, uncert);
[data_bal] = RecallOdeFcn(@descent_ballistic, Td, Yd, settings, uw, vw, ww, uncert);
data_bal.state.Y = Yd;
data_bal.state.T = Td;
save('descent_plot.mat', 'data_bal');

Yf = [Ya(:, 1:13); Yd];
Tf = [Ta; Td];

%% TIME, POSITION AND VELOCITY AT APOGEE
% Usefull values for the plots
bound_value.t = Ta(end);
bound_value.X = [Ya(end, 2), Ya(end, 1), -Ya(end, 3)];
bound_value.V = quatrotate(quatconj(Ya(end, 10:13)), Ya(end, 4:6));

