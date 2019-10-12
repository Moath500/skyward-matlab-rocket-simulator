function [Tf, Yf, Ta, Ya, bound_value] = std_run(settings)
%{ 

STD_RUN - This function runs a standard (non-stochastic) simulation

INTPUTS: 
            - settings, rocket data structure;

OUTPUTS:
            - Tf, Total integration time vector; 
            - Yf, Total State Matrix;
            - Ta, Ascent Integration time vector;
            - Ya, Ascent State Matrix;
            - bound_value, Usefull values for the plots;

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

if settings.para(settings.Npara).z_cut ~= 0 
    error('The landing will be not achived, check the final altitude of the last parachute in config.m')
end

%% WIND GENERATION
if settings.wind.model || settings.wind.input   % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0;
else
    [uw, vw, ww, Azw] = wind_const_generator(settings.wind.AzMin, settings.wind.AzMax,...
        settings.wind.ElMin, settings.wind.ElMax, settings.wind.MagMin, settings.wind.MagMax);
    
    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
    
end

if settings.wind.input && all(settings.wind.input_uncertainty ~= 0)
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
    
    uncert = rand(1, 2).*unc;
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
Q0 = angle2quat(settings.PHI, settings.OMEGA, 0*pi/180, 'ZYX')';

%% ASCENT
% ascent phase computation
Y0a = [X0; V0; W0; Q0; settings.m0; settings.Ixxf; settings.Iyyf; settings.Izzf];
[Ta,Ya] = ode113(@ascent, [0, tf], Y0a, settings.ode.optionsasc1, settings, uw, vw, ww, uncert); % till the apogee

if settings.para(1).delay ~= 0 % checking if the actuation delay is different from zero
    [Ta2,Ya2] = ode113(@ascent, [Ta(end), Ta(end) + settings.para1.delay], Ya(end,:), settings.ode.optionsasc2,...
        settings, uw, vw, ww, uncert); % till end of the delay
    
    Ta = [Ta; Ta2(2:end ) ];
    Ya = [Ya; Ya2(2:end, :)];
end

[data_ascent] = RecallOdeFcn(@ascent, Ta, Ya, settings, uw, vw, ww, uncert);
data_ascent.state.Y = Ya;
data_ascent.state.T = Ta;
save('ascent_plot.mat', 'data_ascent');

%% PARATCHUTES
% Initial Condition are the last from ascent (need to rotate because
% velocities are in body axes)
Y0p = [Ya(end, 1:3) quatrotate(quatconj(Ya(end, 10:13)),Ya(end, 4:6))];
data_para = cell(settings.Npara, 1);
Yf = Ya(:, 1:6);
Tf  = Ta;
t0p = Ta(end);

for i = 1:settings.Npara
    para = i;
    [Tp, Yp] = ode113(@descent_parachute, [t0p, tf], Y0p, settings.ode.optionspara,...
        settings, uw, vw, ww, para, uncert);

    [data_para{para}] = RecallOdeFcn(@descent_parachute, Tp, Yp, settings, uw, vw, ww, para, uncert);
    data_para{para}.state.Y = Yp;
    data_para{para}.state.T = Tp;
    
    % total state
    Yf = [Yf; Yp];
    Tf = [Tf; Tp];
    
    % updating ODE starting conditions
    Y0p = Yp(end, :);
    t0p = Tp(end);
    
end
    
save('descent_para_plot.mat', 'data_para')
    
%% TIME, POSITION AND VELOCITY AT PARACHUTES DEPLOYMENT
% Usefull values for the plots
bound_value = struct;
bound_value(1).t = Ta(end);
bound_value(1).X = [Ya(end, 2), Ya(end, 1), -Ya(end, 3)];
bound_value(1).V = quatrotate(quatconj(Ya(end, 10:13)), Ya(end, 4:6));

for i = 1:settings.Npara
    bound_value(i+1).t = data_para{i}.state.T(end);
    bound_value(i+1).X = [data_para{i}.state.Y(end, 2), data_para{i}.state.Y(end, 1), -data_para{i}.state.Y(end, 3)];
    bound_value(i+1).V = [data_para{i}.state.Y(end, 4), data_para{i}.state.Y(end, 5), -data_para{i}.state.Y(end, 6)];
end



