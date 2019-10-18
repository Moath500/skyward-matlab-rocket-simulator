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

% launchpad 
settings.z0 = 1416;                   %[m] Launchpad Altitude
settings.lrampa = 5.3;                %[m] LaunchPad route (launchpad length-distance from ground of the first hook)

% launchpad directions
% for a single run the maximum and the minimum value of the following
% angles must be the same.
settings.OMEGA = 85*pi/180;        %[rad] Minimum Elevation Angle, user input in degrees (ex. 80)       
settings.PHI = 0*pi/180;       %[rad] Maximum Azimuth Angle from North Direction, user input in degrees (ex. 90)

%% ENGINE DETAILS

engine = 2;

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
        settings.motor.Name = 'K550';
        settings.motor.exp_time = [0 0.13 0.38 0.63 0.88 1.14 1.39...
            1.64 1.9 2.15 2.40 2.66 2.91 3.16 3.5];                         %[s]
        
        settings.motor.exp_thrust = [ 0 139.8 158.07 171.978 178.769 ...
            178.247 158.859 132.922 111.005 92.082 74.075 44.837 16.156...
            4.589 0.000  ] * 9.81/2.2;                                      % [N]
        
        settings.ms = 6.330;                                                % [kg]   Structural Mass
        settings.mp = 0.889;                                                % [kg]   Propellant Mass
        settings.m0 = settings.ms + settings.mp;                            % [kg]   Overall Mas                                       
        settings.mnc = 0.300;                                               % [kg]   Nosecone Mass
        settings.tb = 3.5;                                                  % [s]    Burning time
        settings.mfr = settings.mp/settings.tb;                             % [kg/s] Mass Flow Rate
        
end

%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.09;                          % [m]      Caliber (Fuselage Diameter)
settings.S = 0.0064;                        % [m^2]    Cross-sectional Surface
L = 2.02;                                   % [m]      Rocket length

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.008795446;                    % [kg*m^2] Inertia to x-axis
settings.Iyyf = 2.050393979;                    % [kg*m^2] Inertia to y-axis
settings.Izzf = 2.050413838;                    % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.008472446;                    % [kg*m^2] Inertia to x-axis
settings.Iyye = 1.712284592;                    % [kg*m^2] Inertia to y-axis
settings.Izze = 1.712304085;                    % [kg*m^2] Inertia to z-axis



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
    'Events',@event_landing);              %ODE options for drogue

settings.ode.optionsdesc = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@event_landing);                   %ODE options for ballistic descent



%% Random wind model


% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin = 10;                   % [m/s] Minimum Magnitude
settings.wind.MagMax = 10;                   % [m/s] Maximum Magnitude
settings.wind.ElMin = 0*pi/180;             % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax = 0*pi/180;             % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin = (180)*pi/180;         % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax = (180)*pi/180;         % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind aziumt angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West
