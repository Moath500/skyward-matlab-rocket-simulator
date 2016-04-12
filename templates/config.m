%CONFIG SCRIPT - This script sets up all the parameters for the simulation (H1 line)
% All the parameters are stored in the "settings" structure.
%

% TODO: GUI to fill automatically this configuration script

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% License: 2-clause BSD

% clear all
% close all
% clc

% ROCKET NAME
settings.rocket_name = 'R2A';

% LAUNCHPAD %
settings.z0 = 1417; %Launchpad Altitude
settings.lrampa = 7; %LaunchPad route (launchpad lenght-distance from ground of the first hook)


%STARTING ATTITUDE SETUP %
settings.OMEGA = 85*pi/180; %Elevation Angle
settings.PHI = 300*pi/180; %Azimuth Angle from North Direction

% ENGINE DETAILS %
% % Thrust on time profile (IV Order Interpolation)
% settings.T_coeffs = [8034.5, 0, 0, 0, 0]; % constant thrust
% data from the motor details sheet (experimental plot of the Thrust)

%Cesaroni PRO 150 White Thunder
%settings.motor.exp_time = [0 0.05 0.15 0.5 0.6 0.74 0.85 1.15 1.7 2.4 3 4 4.5 4.8 4.9 5 5.05 5.1 5.15 5.2]; %s
%settings.motor.exp_thrust = [0 8605.1 7900 8400 8400 8250 8200 8300 8400 8400 8200 7800 7600 7450 7350 7300 4500 500 100 0]; %N

%Cesaroni PRO 150 SkidMark
settings.motor.exp_time = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.2 1.8 3.2 3.6 4.8 6 7 7.2 7.6 7.8 7.9 8 8.1 8.19]; %s
settings.motor.exp_thrust = [0 3400 3100 3000 3300 3400 3500 3700 3700 3800 4000 4081.6 3900 3800 3700 3500 3350 3200 3000 2000 750 0]; %N


settings.m0 = 60.29; %kg Overall Mass
settings.ms = 43.13; %kg Structural Mass
settings.mp = settings.m0 - settings.ms; %kg Propellant Mass
settings.tb = 8.19; %sec Burning Time
settings.mfr = settings.mp/settings.tb; %kg/sec Mass Flow Rate

% GEOMETRY DETAILS %
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.174; %m Caliber (Fuselage Diameter)
settings.S = 0.02378; %m2 Cross-sectional Surface
L = 4.396; %m lunghezza missile

% MASS GEOMERTY DETAILS %
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% % NB: le inerzie vengono inserite nel codice di "apogee_reached.m"
% % HP: inerzie razzo = inerzie cilindro pieno circoscritto
% settings.Ixxf=settings.m0*(settings.C/2)^2/2; %Inertia to x-axis (Full)
% settings.Ixxe=settings.ms*(settings.C/2)^2/2; %Inertia to x-axis (Empty)
% settings.Iyyf=settings.m0.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Full)
% settings.Iyye=settings.ms.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Empty)
% settings.Izzf=settings.Iyyf; %Inertia to z-axis (Full)
% settings.Izze=settings.Iyye; %Inertia to z-axis (Empty)


% prima approssimazione inerzie
settings.Ixxf= 0.22723;  %kg*m2 Inertia to x-axis (Full)
settings.Ixxe= 0.16584;  %kg*m2 Inertia to x-axis (Empty)
settings.Iyyf= 50.01003; %kg*m2 Inertia to y-axis (Full)
settings.Iyye= 37.51208; %kg*m2 Inertia to y-axis (Empty)
settings.Izzf= 50.01116; %kg*m2 Inertia to z-axis (Full)
settings.Izze= 37.51321; %kg*m2 Inertia to z-axis (Empty)

% AERODYNAMICS DETAILS %
% This coefficients are intended to be obtained through MISSILE DATCOM
% (than parsed with the tool datcom_parser.py)
% The files are stored in the ../data folder with a naming convention of 
% rocket_name_full.mat | rocket_name_empty.mat
% e.g. R1X_full.mat etc..

%Relative Path of the data files (default: ../data/). Remember the trailing
% slash!!

DATA_PATH = '../data/'; 
filename = strcat(DATA_PATH, settings.rocket_name);

%Coefficients in full configuration (with all the propellant embarqued)
filename_full = strcat(filename,'_full.mat');
CoeffsF = load(filename_full,'Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');

%Coefficients in empty configuration (all the propellant consumed)
filename_empty = strcat(filename,'_empty.mat');
CoeffsE = load(filename_empty,'Coeffs');
settings.CoeffsE = CoeffsE.Coeffs;
clear('CoeffsE');

%Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
%empty and full configuration
s = load(filename_full,'State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
clear('s');

%PARACHUTES DETAILS %
%%% DROGUE 1 %%%
settings.para1.S = 1.153; %m2 %Surface
settings.para1.mass = 0.0692; %kg %Parachute Mass
settings.para1.CD = 0.9; %Parachute Drag Coefficient
settings.para1.CL = 0; %Parachute Lift Coefficient

%Altitude of Drogue 2 Opening
settings.zdrg2 = 2000;

%%% DROGUE 2 %%%
settings.para2.S = 11.520; %m2 %Surface
settings.para2.mass = 0.691; %kg %Parachute Mass
settings.para2.CD = 0.59; %Parachute Drag Coefficient
settings.para2.CL = 0; %Parachute Lift Coefficient

%Altitude of Main Parachute Opening
settings.zmain = 1000;

%%% MAIN - ROGALL %%%
%The drogue parachute effects are neglected

settings.para3.S = 6.327; %m2 %Surface
settings.para3.mass = 0.380; %kg %Parachute Mass
settings.para3.CD = 0.4; %Parachute Mass
settings.para3.CL = 0.9; %Parachute Lift Coefficient

% INTEGRATION OPTIONS %
settings.ode.timeasc = 0:0.01:2000; %Time span for ascend 
settings.ode.timedrg1 = 0:0.01:2000; %Time span for drogue 1
settings.ode.timedrg2 = 0:0.01:2000; %Time span for drogue 2
settings.ode.timemain = 0:0.01:2000; %Time span for main (rogall)

settings.ode.optionsasc = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@apogee,'InitialStep',1); %ODE options for ascend

settings.ode.optionsdrg1 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@drg2_opening); %ODE options for drogue

settings.ode.optionsdrg2 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@main_opening); %ODE options for drogue

settings.ode.optionsmain = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@crash); %ODE options for ascend


% WIND DETAILS %
% Wind is randomly generated. Setting the same values for min and max will
% fix the parameters of the wind.

settings.wind.MagMin = 1.646; %Minimum Magnitude
settings.wind.MagMax = 1.646; %Maximum Magnitude
settings.wind.ElMin = 0*pi/180; %Minimum Elevation
settings.wind.ElMax = 0*pi/180; %Maximum Elevation (Max == 90 Deg)
settings.wind.AzMin = (180 - (39))*pi/180; %Minimum Azimuth
settings.wind.AzMax = (180 - (39))*pi/180; %Maximum Azimuth

% NOTA: angolo di azimut vento = 180 + ...
% significa che ho il vento contro con un angolo di ...
% 180+(..) diretto verso sud-ovest
% 180-(..) diretto verso sud-est

settings.wind.Lat = 41.846647; %Latitude of launching site (Roccaraso)
settings.wind.Long = 14.078528; %Longitude of launching site
settings.wind.Day = randi(365); %Day of the launch (random)
settings.wind.Seconds = randi(86400); %Second of the day (random)

% settings.wind.ww = vert_windgen(settings.wind.MagMin,settings.wind.MagMax); %Vertical wind speed
settings.wind.ww = 0; % no vertical wind

% STOCHASTIC DETAILS %
%If N>1 the stochastic routine is fired (different standard plots)
settings.stoch.N = 1; %Number of iterations
settings.stoch.parallel = false; %Using parallel or not parallel

% PLOT DETAILS %
settings.plot = true; % Set to True to Plot with default plots
settings.tSteps = 250; % Set the number of time steps to visualize
settings.DefaultFontSize = 10; %Default font size for plot
settings.DefaultLineWidth = 1; % Default Line Width for plot