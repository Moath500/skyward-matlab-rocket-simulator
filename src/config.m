% CONFIG - This script sets up all the parameters for the simulation 
% All the parameters are stored in the "settings" structure.

%% LAUNCH SETUP

% rocket name
settings.rocket_name = "R2A_hermes_V1";

% launchpad 
settings.z0 = 1416;                   %[m] Launchpad Altitude
settings.lrampa = 5.3;                %[m] LaunchPad route (launchpad length-distance from ground of the first hook)
settings.lat0 = 41.809017;                                                          % Launchpad latitude
settings.lon0 = 14.054264;                                                          % Launchpad longitude
settings.funZ = funZ_gen('zdata.mat',settings.lat0,settings.lon0,true,'xy');        % Altitude map computation

% launchpad directions
% for a single run the maximum and the minimum value of the following
% angles must be the same.
settings.OMEGAmin = 85*pi/180;        %[rad] Minimum Elevation Angle, user input in degrees (ex. 80)
settings.OMEGAmax = 85*pi/180;        %[rad] Maximum Elevation Angle, user input in degrees (ex. 80)
settings.PHImin = 0*pi/180;           %[rad] Minimum Azimuth Angle from North Direction, user input in degrees (ex. 90)
settings.PHImax = 0*pi/180;           %[rad] Maximum Azimuth Angle from North Direction, user input in degrees (ex. 90)
settings.upwind = false;              % If true, phi is selected according to wind direction (constant wind model only)
settings.PHIsigma = 0*pi/180;         % Stocasthic simulation only

%% ENGINE DETAILS

% Aerotech K550W-L

settings.motor.Name = 'K550';
settings.motor.exp_time = [0 0.13 0.38 0.63 0.88 1.14 1.39...
    1.64 1.9 2.15 2.40 2.66 2.91 3.16 3.5];                         %[s]

settings.motor.exp_thrust = [ 0 139.8 158.07 171.978 178.769 ...
    178.247 158.859 132.922 111.005 92.082 74.075 44.837 16.156...
    4.589 0.000  ] * 9.81/2.2;                                      % [N]

settings.mp = 0.889;                                                % [kg]   Propellant Mass
settings.mnc = 0.300;                                               % [kg]   Nosecone Mass
settings.tb = 3.5;                                                  % [s]    Burning time
settings.mfr = settings.mp/settings.tb;                             % [kg/s] Mass Flow Rate
settings.ms = 4.8;                                                  % [kg]   Total Mass
settings.m0 = settings.ms + settings.mp;                            % [kg]   Structural Mass

%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.09;                          % [m]      Caliber (Fuselage Diameter)
settings.S = 0.0064;                        % [m^2]    Cross-sectional Surface
L = 1.895;                                  % [m]      Rocket length

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.004612382;                    % [kg*m^2] Inertia to x-axis
settings.Iyyf = 1.194050037;                    % [kg*m^2] Inertia to y-axis
settings.Izzf = 1.194116615;                    % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.004293335;                    % [kg*m^2] Inertia to x-axis
settings.Iyye = 0.931311998;                    % [kg*m^2] Inertia to y-axis
settings.Izze = 0.931376985;                    % [kg*m^2] Inertia to z-axis
     
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
settings.para(1).S = 0.7;              % [m^2]   Surface
settings.para(1).mass = 0.075;         % [kg]   Parachute Mass
settings.para(1).CD = 0.75;            % [/] Parachute Drag Coefficient
settings.para(1).CL = 0;               % [/] Parachute Lift Coefficient
settings.para(1).delay = 0;            % [s] drogue opening delay 

% drogue 2
settings.para(2).S = 10.3;             % [m^2]   Surface
settings.para(2).mass = 0.45;          % [kg]   Parachute Mass
settings.para(2).CD = 0.7;             % [/] Parachute Drag Coefficient
settings.para(2).CL = 0.0;             % [/] Parachute Lift Coefficient
settings.para(2).z_opening = 300;      % [m] Altitude of drogue opening


%% INTEGRATION OPTIONS

settings.ode.final_time =  2000;                 % [s] Final integration time

% create an option structure for the integrations:

% - AbsTol is the threshold below which the value of the solution becomes unimportant
% - RelTol is the tolerance betweeen two consecutive values
% - Events is the event function that defines when the integration must be
% - stopped (it has to be created)
% - InitialStep is the highest value tried by the solver

settings.ode.optionsasc1 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_apogee,'InitialStep',1);    %ODE options for ascend

settings.ode.optionsasc2 = odeset('AbsTol',1E-3,'RelTol',1E-3,'InitialStep',1);    
%ODE options for balistic descent between the apogee and the first drogue opening 

settings.ode.optionspara = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_para_opening);              %ODE options for drogue

settings.ode.optionsrog = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_landing);              %ODE options for descent

settings.ode.optionsdesc = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@event_landing);                   %ODE options for ballistic descent


%% WIND DETAILS

% select which model you want to use:

%%%%% Matlab Wind Model
settings.wind.model = false;
% matlab hswm model, wind model on altitude based on historical data

% input Day and Hour as arrays to run stochastic simulations

settings.wind.DayMin = 105;                         % [d] Minimum Day of the launch 
settings.wind.DayMax = 105;                         % [d] Maximum Day of the launch
settings.wind.HourMin = 13;                         % [h] Minimum Hour of the day
settings.wind.HourMax = 13;                         % [h] Maximum Hour of the day
settings.wind.ww = 0;                               % [m/s] Vertical wind speed
                

%%%%% Input wind 
settings.wind.input = false;
% Wind is generated for every altitude interpolating with the coefficient defined below

% first row: wind magnitude [m/s]
% secon row: wind azimut angle (toward wind incoming direction) [deg]
% third row: altitude

% POST INTEGRATION
% settings.wind.input_matr = [ 9*ones(1,7)
%     337.5*ones(1,7)
%     0    100  600  750  900  1500 2500 ];

% IN RAMPA
V0 = 3;
C = [0 0 10 15 20 30 40];
    
settings.wind.input_matr = [ (V0+V0*C/100)
    120*ones(1,7)
    0    100  600  750  900  1500 2500 ];

settings.wind.input_uncertainty = [30,20];
% settings.wind.input_uncertainty = [a,b];      wind uncertanties:
% - a, wind magnitude percentage uncertanty: magn = magn *(1 +- a)
% - b, wind direction band uncertanty: dir = dir 1 +- b

%%%%% Random wind model
% if both the model above are false

% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin = 4;                 % [m/s] Minimum Magnitude
settings.wind.MagMax = 4;                  % [m/s] Maximum Magnitude
settings.wind.ElMin = 0*pi/180;             % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax = 0*pi/180;             % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin = (270)*pi/180;           % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax = (270)*pi/180;         % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind aziumt angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West

%% BALLISTIC SIMULATION
% Set to True to run a ballistic (without drogues) simulation

settings.ballistic = false;    

%% STOCHASTIC DETAILS
% If N > 1 the stochastic routine is started

settings.stoch.N = 1;                             % Number of cases

%%% launch probability details
settings.stoch.prob.x_lim = 2e3;                  % Max ovest displacement [m]
settings.stoch.prob.V_lim = 50;                   % Max drogue velocity [Pa]

%%% Safe Ellipse
settings.prob.SafeEllipse.a = 1100;
settings.prob.SafeEllipse.b = 2800;
settings.prob.SafeEllipse.x0  = 0;
settings.prob.SafeEllipse.y0 = -300;
settings.prob.SafeEllipse.alpha = 10;

%% PLOT DETAILS

settings.plots = true;
settings.terrain = true;

%% LANDING POINTS
settings.landing_map = true;
settings.map_file = 'map_roccaraso.jpg'; % name of map for landing points
settings.map_xaxis = [-5000 5000];  % limits for the data of the landing map
settings.map_yaxis = [-5000 5000];