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
settings.rocket_name = "R2A_hermes";

% launchpad 6
settings.z0 = 1416;                %[m] Launchpad Altitude
settings.lrampa = 4;               %[m] LaunchPad route (launchpad length-distance from ground of the first hook)


% starting altitude
settings.OMEGA = 80*pi/180;        %[rad] Elevation Angle, user input in degrees (ex. 80)
settings.PHI = 0*pi/180;           %[rad] Azimuth Angle from North Direction, user input in degrees (ex. 90)

%% ENGINE DETAILS

% sintax:
% engine 1 -> Aerotech K1000-T
% engine 2 -> Aerotech K560-W
engine = 1;

switch engine
       case 1 
        settings.motor.Name = 'K1000T-P';
        settings.motor.exp_time =   [0   0.0150    0.0250    0.0950  ...
            0.2000    0.3000    0.4000    0.5000    0.6000    0.7000 ...
            0.8000    0.9000    1.0000    1.1000    1.2000    1.3000 ...
            1.4000    1.5000    1.6000    1.7000    1.8000    1.9000 ...
            2.0000    2.1000    2.1800    2.2000    2.2180    2.2690 ...  
            2.3000    2.3320    2.3560    2.3890    2.4360    2.5000];      % [s]
        
        settings.motor.exp_thrust = [ 0.8951    1.1198    1.0933    1.0966...
            1.1099    1.1165    1.1231    1.1330    1.1396    1.1363 ...
            1.1363    1.1363    1.1396    1.1330    1.1297    1.1264 ...
            1.1198    1.1099    1.0966    1.0636    1.0174    0.9711 ...
            0.9150    0.8687    0.8654    0.8786    0.8588    0.6705 ...
            0.5780    0.4459    0.3369    0.2246    0.1057     0].*1e3  ;   % [N]
        
        settings.m0 = 8.545;                                                % [kg]   Overall Mass
        settings.ms = 7.177;                                                % [kg]   Structural Mass
        settings.mnc = 0.120;                                               % [kg]   Nosecone Mass
        settings.tb = 2.5;                                                  % [s]    Burning time
        settings.mp = settings.m0-settings.ms;                              % [kg]   Propellant Mass
        settings.mfr = settings.mp/settings.tb;                             % [kg/s] Mass Flow Rate
   case 2
        settings.motor.Name = 'K560W-PS';
        settings.motor.exp_time = [ 0   0.2900    0.4840    0.6790 ...
            0.8740    1.0690    1.2640    1.4580  1.6530    1.8480 ...
            2.0430    2.2370    2.4320    2.6270  2.8220    3.0170 ...
            3.2110    3.4060    3.6010    3.7960  3.9910    4.1850 ... 
            4.3800    4.5750    4.7700    4.9650];                          %[s]
        
        setting.motor.exp_thrust = [552.1230  645.4030  681.1090 ...
               716.1670  742.6780  764.7780  775.7100  785.8590  ...
               789.3050  789.0770  744.6220  676.8860  614.7110  ...
               557.9080  503.6410  455.5040  412.0450  372.9630  ...
               335.9870  307.3460  279.8560  223.4910   70.4410  ...
               10.0280   2.4450       0];                                   % [N]
        
        settings.m0 = 8.558;                                                % [kg]   Overall Mas
        settings.ms = 7.134;                                                % [kg]   Structural Mass
        settings.mp = settings.m0-settings.ms;                              % [kg]   Propellant Mass
        settings.mnc = 0.149;                                               % [kg]   Nosecone Mass
        settings.tb = 4.9;                                                  % [s]    Burning time
        settings.mfr = settings.mp/settings.tb;                             % [kg/s] Mass Flow Rate
end

%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.09;                          % [m]      Caliber (Fuselage Diameter)
settings.S = 0.0064;                        % [m^2]    Cross-sectional Surface
L = 1.91;                                   % [m]      Rocket length

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.01101;                    % [kg*m^2] Inertia to x-axis
settings.Iyyf = 2.1137;                     % [kg*m^2] Inertia to y-axis
settings.Izzf = 2.1139;                     % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.01024;                    % [kg*m^2] Inertia to x-axis
settings.Iyye = 1.76211;                    % [kg*m^2] Inertia to y-axis
settings.Izze = 1.76235;                    % [kg*m^2] Inertia to z-axis


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
settings.para1.S = 1.55;             % [m^2]   Surface
settings.para1.mass = 0.25;          % [kg]   Parachute Mass
settings.para1.CD = 0.8;             % [/] Parachute Drag Coefficient
settings.para1.CL = 0;               % [/] Parachute Lift Coefficient

% rogallo wing
settings.para2.S = 17.5;             % [m^2]   Surface
settings.para2.mass = 1.140;         % [kg]   Parachute Mass
settings.para2.CD = 0.59;            % [/] Parachute Drag Coefficient
settings.para2.CL = 0;               % [/] Parachute Lift Coefficient
settings.zdrg2 = 1000;               % [m] Altitude of drogue 2 opening


%% INTEGRATION OPTIONS

settings.ode.timeasc =  [0 100];     % [s]   Time span for ascend
settings.ode.timedrg1 = [0 400];     % [s]   Time span for drogue 1
settings.ode.timedrg2 = [0 400];     % [s]   Time span for drogue 2
settings.ode.timedesc = [0 100];     % [s]   Time span for ballistic descent



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
    'Events',@event_landing);              %ODE options for drogue

settings.ode.optionsdesc = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@event_landing);                   %ODE options for ballistic descent


%% WIND DETAILS


% Settings for the Wind Model

settings.wind.model = false;
% set to true for hwsm wind model
% set to false for random wind model
% input Day and Hour as arrays to run stochastic simulations

settings.wind.Lat = 41.800833;                      % [deg] Latitude of launching site
settings.wind.Long = 14.055833;                      % [deg] Longitude of launching site
settings.wind.DayMin = 105;                         % [d] Minimum Day of the launch 
settings.wind.DayMax = 105;                         % [d] Maximum Day of the launch
settings.wind.HourMin = 13;                          % [h] Minimum Hour of the day
settings.wind.HourMax = 13;                         % [h] Maximum Hour of the day
settings.wind.ww = 0;                               % [m/s] Vertical wind speed
                

% Input wind 

% Wind is generated for every altitude interpolating with the coefficient defined below

% first row: wind magnitude [m/s]
% secon row: wind azimut angle [deg]
% third row: altitude

settings.wind.input = false;
settings.wind.input_matr = [ 5    7    9   10    11    11   13   12   13  
                             250  260  260 260   260   260  270  270  270   
                             0    100  600 750   900   1500 2000 3000 6000 ];
                         
settings.wind.input_uncertainty = 0;              % [perc] uncertainty percentage

% Random wind model

% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin = 5;                   % [m/s] Minimum Magnitude
settings.wind.MagMax = 5;                   % [m/s] Maximum Magnitude
settings.wind.ElMin = 0*pi/180;              % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax = 0*pi/180;              % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin = (90)*pi/180;           % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax = (90)*pi/180;           % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind aziumt angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West

%% BALLISTIC SIMULATION

settings.ballistic = true;                  % Set to True to run a standard ballistic (without drogues) simulation

%% SECOND DROGUE FAILURE
% simulation in which second drogue does not open and thus landing is
% achieved in ballistic (settings.ballistic must be true)

settings.sdf = false;

%% APOGEE ONLY
% simulation stopped when reaching the apogee and thus there is no
% descend phase.   Only available for standard stochastic runs !!!

settings.ao = false;

%% STOCHASTIC DETAILS
% If N > 1 the stochastic routine is started

settings.stoch.N = 1;             % Number of cases

%% PLOT DETAILS

settings.plots = true;
settings.only_XCP = false; % plot only the stability margin
settings.landing_map = false;
