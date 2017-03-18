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

% LAUNCHPAD %
settings.z0 = 1200; %Launchpad Altitude
settings.lrampa = 3.5; %LaunchPad Length


%STARTING ATTITUDE SETUP %
settings.OMEGA = 89*pi/180; %Elevation Angle
settings.PHI = 0*pi/180; %Azimuth Angle from North Direction

% ENGINE DETAILS %
% settings.T_coeffs = fliplr([-332.799894219468, 2209.93106349784, -5265.48638128018, 4608.82884815720, -133.133036494531]);
                    %[0 438.918868 -185.329879 33.642643 -2.517265]; 
%Thrust on time profile (IV Order Interpolation)
settings.T_coeffs = [8034.5, 0, 0, 0, 0];

settings.ms = 40.4; %kg Dry Mass (No Propellant)
settings.tb = 5.12; %sec Burning Time
settings.mp = 18.6; %kg Mass of Propellant
settings.mfr = settings.mp/settings.tb; %kg/sec Mass Flow Rate
settings.m0 = settings.ms+settings.mp; %kg Total Starting Mass

% GEOMETRY DETAILS %
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.174; %m Caliber (Fuselage Diameter)
settings.S = 0.02378; %m2 Cross-sectional Surface
L = 4.250; %m lunghezza missile

% MASS GEOMERTY DETAILS %
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% NB: le inerzie vengono inserite nel codice di "apogee_reached.m"
% HP: inerzie razzo = inerzie cilindro pieno circoscritto
settings.Ixxf=settings.m0*(settings.C/2)^2/2; %Inertia to x-axis (Full)
settings.Ixxe=settings.ms*(settings.C/2)^2/2; %Inertia to x-axis (Empty)
settings.Iyyf=settings.m0.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Full)
settings.Iyye=settings.ms.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Empty)
settings.Izzf=settings.Iyyf; %Inertia to z-axis (Full)
settings.Izze=settings.Iyye; %Inertia to z-axis (Empty)

% % prova senza inerzie
% settings.Ixxf=0.21; %Inertia to x-axis (Full)
% settings.Ixxe=0.14; %Inertia to x-axis (Empty)
% settings.Iyyf=20; %Inertia to y-axis (Full)
% settings.Iyye=15; %Inertia to y-axis (Empty)
% settings.Izzf=20; %Inertia to z-axis (Full)
% settings.Izze=15; %Inertia to z-axis (Empty)

% AERODYNAMICS DETAILS %
% This coefficients are intended to be obtained through MISSILE DATCOM
% (than parsed with the tool datcom_parser.py)
CoeffsF = load('for006_full.mat','Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');
CoeffsE = load('for006_empty.mat','Coeffs');
settings.CoeffsE = CoeffsE.Coeffs;
clear('CoeffsE');

%Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
%empty and full configuration
s = load('for006_full.mat','State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
clear('s');

%PARACHUTES DETAILS %
%%% DROGUE %%%
settings.para1.S = 0.52; %Surface
settings.para1.mass = 0.02; %Parachute Mass
settings.para1.CD = 0.9; %Parachute Drag Coefficient


%Altitude of Main Parachute Opening
settings.zmain = 250;

%%% MAIN %%%
%The drogue parachute effects are neglected

settings.para2.S = 11.65; %Surface
settings.para2.mass = 0.43; %Parachute Mass
settings.para2.CD = 0.59; %Parachute Mass


% INTEGRATION OPTIONS %
settings.ode.timeasc = 0:0.01:2000; %Time span for ascend 
settings.ode.timedrg = 0:0.01:2000; %Time span for drogue 
settings.ode.timemain = 0:0.01:2000; %Time span for main 

settings.ode.optionsasc = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@apogee,'InitialStep',1); %ODE options for ascend
settings.ode.optionsdrg = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@main_opening); %ODE options for drogue
settings.ode.optionsmain = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@crash); %ODE options for ascend

% WIND DETAILS %
% Wind is randomly generated. Setting the same values for min and max will
% fix the parameters of the wind.

settings.wind.MagMin = 10; %Minimum Magnitude
settings.wind.MagMax = 10; %Maximum Magnitude
settings.wind.ElMin = 0*pi/180; %Minimum Elevation
settings.wind.ElMax = 0*pi/180; %Maximum Elevation (Max == 90 Deg)
settings.wind.AzMin = 45*pi/180; %Minimum Azimuth
settings.wind.AzMax = 45*pi/180; %Maximum Azimuth

% STOCHASTIC DETAILS %
%If N>1 the stochastic routine is fired (different standard plots)
settings.stoch.N = 1; %Number of iterations
settings.stoch.parallel = true; %Using parallel or not parallel

% PLOT DETAILS %
settings.plot = false; % Set to True to Plot with default plots
settings.tSteps = 250; % Set the number of time steps to visualize
settings.DefaultFontSize = 10; %Default font size for plot
settings.DefaultLineWidth = 1; % Default Line Width for plot