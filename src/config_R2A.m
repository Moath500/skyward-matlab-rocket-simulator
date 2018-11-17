% CONFIG - This script sets up all the parameters for the simulation 
% All the parameters are stored in the "settings" structure.

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% License: 2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

%% LAUNCH SETUP

% rocket name
settings.rocket_name = "R2A";

% launchpad 6
settings.z0 = 5;                                 %[m] Launchpad Altitude
settings.lrampa = 5.5;                           %[m] LaunchPad route (launchpad length-distance from ground of the first hook)


% starting altitude
settings.OMEGA = 80*pi/180;                      %[rad] Elevation Angle, user input in degrees (ex. 80)
settings.PHI = 90*pi/180;                        %[rad] Azimuth Angle from North Direction, user input in degrees (ex. 90)

%% ENGINE DETAILS

% sintax:
% engine = 1 -> Cesaroni PRO 150 White Thunder
% engine = 2 -> Cesaroni PRO 150 SkidMark
% engine = 3 -> Cesaroni PRO 150 BlueStreak
engine = 3;

switch engine
    case 1
        % Cesaroni PRO 150 White Thunder
        % Sampling for thrust interpolation
        settings.motor.Name = 'Cesaroni PRO 150 White Thunder';
        settings.motor.exp_time = [0 0.05 0.15 0.5 0.6 0.74 0.85 1.15 1.7 2.4 3 ...
            4 4.5 4.8 4.9 5 5.05 5.1 5.15 5.2];  % [s]
        settings.motor.exp_thrust = [8605.1 8900 7900 8400 8400 8250 8200 8300 ...
            8400 8400 8200 7800 7600 7450 7350 7300 4500 500 100 0]; % [N]
        
        settings.m0 = 67.761;                    % [kg]    Overall Mass (Burnout + Propellant)
        settings.ms = 43.961;                    % [kg]    Structural Mass (Burnout - Nosecone)
        settings.mp = 18.6;                      % [kg]    Propellant Mass
        settings.tb = 5.12;                      % [s]     Burning Time
        settings.mfr = settings.mp/settings.tb;  % [kg/s]  Mass Flow Rate
    case 2
        % Cesaroni PRO 150 SkidMark
        % Sampling for thrust interpolation
        settings.motor.Name = 'Cesaroni PRO 150 SkidMark';
        settings.motor.exp_time = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.2 1.8 3.2 ...
            3.6 4.8 6 7 7.2 7.6 7.8 7.9 8 8.1 8.19]; % [s]
        settings.motor.exp_thrust = [0 3400 3100 3000 3300 3400 3500 3700 3700 ...
            3800 4000 4081.6 3900 3800 3700 3500 3350 3200 3000 2000 750 0]; % [N]
        
        settings.m0 = 64.9;                      % [kg]    Overall Mass
        settings.ms = 46.8;                      % [kg]    Structural Mass (Burnout)
        settings.mp = settings.m0-settings.ms;   % [kg]    Propellant Mass
        settings.tb = 8.19;                      % [s]     Burning Time
        settings.mfr = settings.mp/settings.tb;  % [kg/s]  Mass Flow Rate
    case 3
        % Cesaroni PRO 150 BlueStreak
        % Sampling for thrust interpolation
        settings.motor.Name = 'Cesaroni PRO 150 BlueStreak';
        settings.motor.exp_time =   [0 0.06 0.1 0.15 0.25 0.45 0.8  1     2    3 ...
            4     5   6   6.8  7.05 7.3 7.6 7.8]; % [s]
        settings.motor.exp_thrust = [0 800 4000 5500 5160 5130 5400 5300 5450 5347 ...
            5160 4950 4700 4400 4400 3800 300 0]; % [N]
        
        settings.m0 = 63.30421;                  % [kg]   Overall Mass
        settings.ms = 44.40621;                  % [kg]   Structural Mass (Burnout)
        settings.mp = settings.m0-settings.ms;   % [kg]   Propellant Mass
        settings.mnc = 6.21;                     % [kg]   Nosecone Mass
        settings.tb = 7.60;                      % [s]    Burning Time
        settings.mfr = settings.mp/settings.tb;  % [kg/s] Mass Flow Rate
end


%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.174;                              % [m]      Caliber (Fuselage Diameter)
settings.S = 0.02378;                            % [m^2]    Cross-sectional Surface
L = 4.43619;                                     % [m]      Rocket length

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.27;                            % [kg*m^2] Inertia to x-axis
settings.Iyyf = 84.42;                           % [kg*m^2] Inertia to y-axis
settings.Izzf = 84.42;                           % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.21;                            % [kg*m^2] Inertia to x-axis
settings.Iyye = 65.77;                           % [kg*m^2] Inertia to y-axis
settings.Izze = 65.77;                           % [kg*m^2] Inertia to z-axis

%% AERODYNAMICS DETAILS
% These coefficients are obtained using MISSILE DATCOM
% (after parsing with the tool datcom_parser.py)
% The files are stored in the ../data folder with the following rule:
% rocket_name_full.mat | rocket_name_empty.mat
% e.g. R2a_full.mat    | R2a_empty.mat
% Relative Path of the data files (default: ../data/). Remember the trailing slash!!

% Coeffs is a 4D matrix given by Datcom that contains the aerodynamics
% coefficient computed for the input parameters (AoA,Betas,Altitudes,Machs)
% Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
% empty and full configuration

DATA_PATH = '../data/';
filename = strcat(DATA_PATH, settings.rocket_name);

% Coefficients in full configuration
filename_full = strcat(filename,'_full.mat');
CoeffsF = load(filename_full,'Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');

% Coefficients in empty configuration
filename_empty = strcat(filename,'_empty.mat');
CoeffsE = load(filename_empty,'Coeffs');
settings.CoeffsE = CoeffsE.Coeffs;
clear('CoeffsE');


s = load(filename_full,'State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
clear('s');


%% PARACHUTES DETAILS

% drogue 1
settings.para1.S = 1.55;                         % [m^2]   Surface
settings.para1.mass = 0.25;                      % [kg]   Parachute Mass
settings.para1.CD = 0.8;                         % [/] Parachute Drag Coefficient
settings.para1.CL = 0;                           % [/] Parachute Lift Coefficient

% drogue 2
settings.para2.S = 17.5;                         % [m^2]   Surface
settings.para2.mass = 1.140;                     % [kg]   Parachute Mass
settings.para2.CD = 0.59;                        % [/] Parachute Drag Coefficient
settings.para2.CL = 0;                           % [/] Parachute Lift Coefficient
settings.zdrg2 = 5000;                           % [m] Altitude of drogue 2 opening

% rogallo wing
% The drogue parachute effects are neglected
settings.para3.S = 15;                           % [m^2]   Surface
settings.para3.mass = 1.466;                     % [kg]   Parachute Mass
settings.para3.CD = 0.4;                         % [/] Parachute Drag Coeff
settings.para3.CL = 0.8;                         % [/] Parachute Lift Coefficient
settings.zrog = 2000;                            % [m] Altitude of Rogallo Opening

%% INTEGRATION OPTIONS

settings.ode.final_time =  2000;                 % [s] Final integration time

% create an option structure for the integrations:

% - AbsTol is the threshold below which the value of the solution becomes unimportant
% - RelTol is the tolerance betweeen two consecutive values
% - Events is the event function that defines when the integration must be
% - stopped (it has to be created)
% - InitialStep is the highest value tried by the solver

settings.ode.optionsasc = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_apogee,'InitialStep',1);    %ODE options for ascend

settings.ode.optionsdrg1 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_drg2_opening);              %ODE options for drogue

settings.ode.optionsdrg2 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_rog_opening);              %ODE options for drogue

settings.ode.optionsrog = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_landing);                   %ODE options for descent

settings.ode.optionsdesc = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_landing);                   %ODE options for ballistic descent


%% WIND DETAILS


% Settings for the Wind Model

settings.wind.model = false;
% set to true for hwsm wind model
% set to false for random wind model
% input Day and Hour as arrays to run stochastic simulations

settings.wind.Lat = 39.552709;                   % [deg] Latitude of launching site
settings.wind.Long = 9.652400;                   % [deg] Longitude of launching site
settings.wind.DayMin = 270;                      % [d] Minimum Day of the launch 
settings.wind.DayMax = 270;                      % [d] Maximum Day of the launch
settings.wind.HourMin = 20;                      % [h] Minimum Hour of the day
settings.wind.HourMax = 20;                      % [h] Maximum Hour of the day
settings.wind.ww = 0;                            % [m/s] Vertical wind speed

% Input wind 

% Wind is generated for every altitude interpolating with the coefficient defined below

% first row: wind magnitude [m/s]
% secon row: wind azimut angle [deg]
% third row: altitude

settings.wind.input = true;
settings.wind.input_matr = [ 5    7    9   10   11   11   13   12   13   13   14   12   10    10
                             250  260  260 260  260  260  270  270  270  270  270  270  270   270
                             0    100  600 750  900  1500 2000 3000 4200 5500 7000 9000 10000 18000];
                         
settings.wind.input_uncertainty = 20;             % [perc] uncertainty percentage
                         
                         

% Random wind model

% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin = 6;                         % [m/s] Minimum Magnitude
settings.wind.MagMax = 10;                        % [m/s] Maximum Magnitude
settings.wind.ElMin = 0*pi/180;                   % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax = 0*pi/180;                   % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin = (90)*pi/180;                % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax = (90)*pi/180;                % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind aziumt angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West

%% BALLISTIC SIMULATION

settings.ballistic = true;                      % Set to True to run a standard ballistic (without drogues) simulation

%% LAST DROGUE FAILURE SIMULATION
% simulation in which rogallo wing does not open and thus landing is
% achieved thanks to the 2nd parachute

settings.ldf = false;

%% APOGEE ONLY
% simulation stopped when reaching the apogee and thus there is no
% descend phase.   Only available for standard stochastic runs!!!

settings.ao = false;

%% STOCHASTIC DETAILS
% If N > 1 the stochastic routine is started

settings.stoch.N = 1;                             % Number of cases
settings.stoch.prob = true;                       % Set to true to compute the launch probability
settings.stoch.x_lim = 2e3;                       % Max ovest displacement [m]
settings.stoch.P_lim = 2e3;                       % Max drogue stress [Pa]

%% PLOT DETAILS

settings.plots = false;
settings.only_XCP = false;                        % plot only the stability margin
settings.landing_map = false;

