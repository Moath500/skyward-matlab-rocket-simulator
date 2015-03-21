%CONFIG SCRIPT - This script sets up all the parameters for the simulation (H1 line)
% All the parameters are stored in the "settings" structure.
%

% TODO: GUI to fill automatically this configuration script

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% License: 2-clause BSD

clear all
close all
clc

% LAUNCHPAD %
settings.z0 = 1200; %Launchpad Altitude
settings.lrampa = 3.5; %LaunchPad Length


%STARTING ATTITUDE SETUP %
settings.OMEGA = 89*pi/180; %Elevation Angle
settings.PHI = 0*pi/180; %Azimuth Angle from North Direction

% ENGINE DETAILS %
settings.T_coeffs =fliplr([-332.799894219468, 2209.93106349784, -5265.48638128018, 4608.82884815720, -133.133036494531]);
                    %[0 438.918868 -185.329879 33.642643 -2.517265]; 
%Thrust on time profile (IV Order Interpolation)
settings.ms = 8; %Dry Mass (No Propellant)
settings.tb = 2.1; %Burning Time
settings.mp = 2; %Mass of Propellant
settings.mfr = settings.mp/settings.tb; %Mass Flow Rate
settings.m0 = settings.ms+settings.mp; %Total Starting Mass

% MASS GEOMERTY DETAILS %
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

settings.Ixxf=0.18; %Inertia to x-axis (Full)
settings.Ixxe=0.17; %Inertia to x-axis (Empty)
settings.Iyyf=3.98; %Inertia to y-axis (Full)
settings.Iyye=3.33; %Inertia to y-axis (Empty)
settings.Izzf=3.98; %Inertia to z-axis (Full)
settings.Izze=3.33; %Inertia to z-axis (Empty)

% GEOMETRY DETAILS %
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 90e-3; %Caliber (Fuselage Diameter)
settings.S = settings.C^2/4*pi; %Cross-sectional Surface

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
settings.Alphas = s.State.Alphas;
settings.Betas = s.State.Betas;
settings.Altitudes = s.State.Altitudes;
settings.Machs = s.State.Machs;
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
settings.ode.optionsdrg = odeset('AbsTol',1E-6,'RelTol',1E-6,...
    'Events',@main_opening); %ODE options for drogue
settings.ode.optionsmain = odeset('AbsTol',1E-6,'RelTol',1E-6,...
    'Events',@crash); %ODE options for ascend

% WIND DETAILS %
% Wind is randomly generated. Setting the same values for min and max will
% fix the parameters of the wind.

settings.wind.MagMin = 5; %Minimum Magnitude
settings.wind.MagMax = 5; %Maximum Magnitude
settings.wind.ElMin = 0*pi/180; %Minimum Elevation
settings.wind.ElMax = 0*pi/180; %Maximum Elevation (Max == 90 Deg)
settings.wind.AzMin = 45*pi/180; %Minimum Azimuth
settings.wind.AzMax = 45*pi/180; %Maximum Azimuth

% STOCHASTIC DETAILS %
%If N>1 the stochastic routine is fired (different standard plots)
settings.stoch.N = 1; %Number of iterations
settings.stoch.parallel = false; %Using parallel or not parallel

% PLOT DETAILS %
settings.plot = true; % Set to True to Plot with default plots
settings.tSteps = 250; % Set the number of time steps to visualize
settings.DefaultFontSize = 20; %Default font size for plot



