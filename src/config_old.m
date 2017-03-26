%CONFIG SCRIPT - This script sets up all the parameters for the simulation (H1 line)
% All the parameters are stored in the "settings" structure.
%

% TODO: GUI to fill automatically this configuration script

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% License: 2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

% clear all
% close all
% clc

% LAUNCHPAD %
settings.z0 = 100;       %Launchpad Altitude
settings.lrampa = 4.443; %LaunchPad route (launchpad lenght-distance from ground of the first hook)


%STARTING ATTITUDE SETUP %
settings.OMEGA = 85*pi/180;     %Elevation Angle
settings.PHI = 270*pi/180;      %Azimuth Angle from North Direction

% ENGINE DETAILS %

% % Thrust on time profile (IV Order Interpolation)
% settings.T_coeffs = [8034.5, 0, 0, 0, 0]; % constant thrust
% data from the motor details sheet (experimental plot of the Thrust)

% SYNTAX:
% engine = 1 -> Cesaroni PRO 150 White Thunder
% engine = 2 -> Cesaroni PRO 150 SkidMark
% engine = 3 -> Cesaroni PRO 150 BlueStreak
engine = 3;
switch engine
    case 1
        % Cesaroni PRO 150 White Thunder
        % Sampling for thrust interpolation
        settings.motor.exp_time = [0 0.05 0.15 0.5 0.6 0.74 0.85 1.15 1.7 2.4 3 ...
            4 4.5 4.8 4.9 5 5.05 5.1 5.15 5.2]; %s
        settings.motor.exp_thrust = [8605.1 8900 7900 8400 8400 8250 8200 8300 ...
            8400 8400 8200 7800 7600 7450 7350 7300 4500 500 100 0]; %N

        settings.m0 = 61.8;                      %kg    Overall Mass
        settings.mp = 18.6;                      %kg    Propellant Mass
        settings.ms = settings.m0 - settings.mp; %kg    Structural Mass
        settings.tb = 5.12;                      %sec   Burning Time
        settings.mfr = settings.mp/settings.tb;  %kg/s  Mass Flow Rate
    case 2
        % Cesaroni PRO 150 SkidMark
        % Sampling for thrust interpolation
        settings.motor.exp_time = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.2 1.8 3.2 ...
            3.6 4.8 6 7 7.2 7.6 7.8 7.9 8 8.1 8.19]; %s
        settings.motor.exp_thrust = [0 3400 3100 3000 3300 3400 3500 3700 3700 ...
            3800 4000 4081.6 3900 3800 3700 3500 3350 3200 3000 2000 750 0]; %N

        settings.m0 = 61.8;                      %kg    Overall Mass
        settings.mp = 17.157;                    %kg    Propellant Mass
        settings.ms = settings.m0 - settings.mp; %kg    Structural Mass
        settings.tb = 8.19;                      %sec   Burning Time
        settings.mfr = settings.mp/settings.tb;  %kg/s  Mass Flow Rate
    case 3
    % Cesaroni PRO 150 BlueStreak
    % Sampling for thrust interpolation
    settings.motor.exp_time =   [0 0.06 0.1 0.15 0.25 0.45 0.8  1     2    3 ...
        4     5   6   6.8  7.05 7.3 7.6 7.8]; %s
    settings.motor.exp_thrust = [0 800 4000 5500 5160 5130 5400 5300 5450 5347 ...
        5160 4950 4700 4400 4400 3800 300 0]; %N
   
    
    settings.mp = 17.7;                      %kg   Propellant Mass
    settings.ms = 39 + 3.9598;               %kg   Structural Mass
    settings.m0 = settings.mp + settings.ms; %kg   Overall Mass without parachutes 
                                                  %considering 3~5 kg DPL compartment
    settings.tb = 7.60;                      %sec  Burning Time
    settings.mfr = settings.mp/settings.tb;  %kg/s Mass Flow Rate
    otherwise
end


% GEOMETRY DETAILS %
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.174;     %m      Caliber (Fuselage Diameter)
settings.S = 0.02378;   %m2     Cross-sectional Surface
L = 4.396;              %m      Rocket length

% MASS GEOMERTY DETAILS %
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% % Note: inerzias are used in "apogee_reached.m"
% % HP: rocket inertias = full finite cilinder inertias
% settings.Ixxf=settings.m0*(settings.C/2)^2/2; %Inertia to x-axis (Full)
% settings.Ixxe=settings.ms*(settings.C/2)^2/2; %Inertia to x-axis (Empty)
% settings.Iyyf=settings.m0.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Full)
% settings.Iyye=settings.ms.*((settings.C/2)^2/4 + L^2/3); %Inertia to y-axis (Empty)
% settings.Izzf=settings.Iyyf; %Inertia to z-axis (Full)
% settings.Izze=settings.Iyye; %Inertia to z-axis (Empty)


% inertias first approximation
settings.Ixxf = 0.14412762;  %kg*m2 Inertia to x-axis (Full)
settings.Ixxe = 0.14412762;  %kg*m2 Inertia to x-axis (Empty)
settings.Iyyf = 30.16891336; %kg*m2 Inertia to y-axis (Full)
settings.Iyye = 30.16891335; %kg*m2 Inertia to y-axis (Empty)
settings.Izzf = 30.16930792; %kg*m2 Inertia to z-axis (Full)
settings.Izze = 30.16930792; %kg*m2 Inertia to z-axis (Empty)

% AERODYNAMICS DETAILS %
% This coefficients are intended to be obtained through MISSILE DATCOM
% (than parsed with the tool datcom_parser.py)
CoeffsF = load('for006_full.mat','Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');
CoeffsE = load('for006_empty.mat','Coeffs');
settings.CoeffsE = CoeffsE.Coeffs;
clear('CoeffsE');

% Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
% empty and full configuration
s = load('for006_full.mat','State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
clear('s');

%PARACHUTES DETAILS %
%%% DROGUE 1 %%%
settings.para1.S = 1.153;           %m2   Surface
settings.para1.mass = 0.0692;       %kg   Parachute Mass
settings.para1.CD = 0.9;            %Parachute Drag Coefficient
settings.para1.CL = 0;              %Parachute Lift Coefficient

%Altitude of Drogue 2 Opening
settings.zdrg2 = 2000;

%%% DROGUE 2 %%%
settings.para2.S = 11.520;          %m2   Surface
settings.para2.mass = 0.691;        %kg   Parachute Mass
settings.para2.CD = 0.59;           %Parachute Drag Coefficient
settings.para2.CL = 0;              %Parachute Lift Coefficient

%Altitude of Main Parachute Opening
settings.zmain = 1000;

%%% MAIN - ROGALLO %%%
%The drogue parachute effects are neglected

settings.para3.S = 6.327;           %m2   Surface
settings.para3.mass = 0.380;        %kg   Parachute Mass
settings.para3.CD = 0.4;            %Parachute Drag Coeff
settings.para3.CL = 0.9;            %Parachute Lift Coefficient

% INTEGRATION OPTIONS %
settings.ode.timeasc = 0:0.01:2000;  %sec   %Time span for ascend 
settings.ode.timedrg1 = 0:0.01:2000; %sec   %Time span for drogue 1
settings.ode.timedrg2 = 0:0.01:2000; %sec   %Time span for drogue 2
settings.ode.timemain = 0:0.01:2000; %sec   %Time span for main (rogallo)
settings.ode.timedesc = 0:0.01:2000; %sec   %Time span for ballistic descent

settings.ode.optionsasc = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@apogee,'InitialStep',1); %ODE options for ascend

settings.ode.optionsdrg1 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@drg2_opening); %ODE options for drogue

settings.ode.optionsdrg2 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@main_opening); %ODE options for drogue

settings.ode.optionsmain = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@crash);        %ODE options for ascend

settings.ode.optionsdesc = odeset('AbsTol',1E-3,'RelTol',1E-12,...
    'Events',@crash);        %ODE options for ballistic descent


% WIND DETAILS %

% Wind is randomly generated. Setting the same values for min and max will
% fix the parameters of the wind.

settings.wind.MagMin = 0;                    %Minimum Magnitude
settings.wind.MagMax = 0;                    %Maximum Magnitude
settings.wind.ElMin = 0*pi/180;              %Minimum Elevation
settings.wind.ElMax = 0*pi/180;              %Maximum Elevation (Max == 90 Deg)
settings.wind.AzMin = (180 + (135))*pi/180;  %Minimum Azimuth
settings.wind.AzMax = (180 + (135))*pi/180;  %Maximum Azimuth

% NOTE: wind aziumt angle = 180 + ...
% means that I'm upwind with an angle of ...
% 180+(..) heading sud-ovest
% 180-(..) heading sud-est

% Settings for the Wind Model
settings.wind.Lat = 52.85;       %Latitude of launching site 
settings.wind.Long = 16.033333;  %Longitude of launching site
settings.wind.Day = 180;         %Day of the launch 
settings.wind.Seconds = 36000;   %Second of the day 

% settings.wind.ww = vert_windgen(settings.wind.MagMin,settings.wind.MagMax); 
%Vertical wind speed
settings.wind.ww = 0; % no vertical wind

% BALLISTIC SIMULATION
settings.ballistic = false;     % Set to True to run a ballistic simulation

% STOCHASTIC DETAILS %
%If N>1 the stochastic routine is fired (different standard plots)
settings.stoch.N = 1;            % Number of iterations
settings.stoch.parallel = false; % Using parallel or not parallel

% PLOT DETAILS %
settings.plot = false;         % Set to True to Plot with default plots
settings.tSteps = 250;         % Set the number of time steps to visualize
settings.DefaultFontSize = 10; % Default font size for plot
settings.DefaultLineWidth = 1; % Default Line Width for plot