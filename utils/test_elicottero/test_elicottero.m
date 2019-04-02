clc; clear; close all


%% CONFIG

optionsdrg1 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_drg2_opening);          
optionsdrg2 = odeset('AbsTol',1E-3,'RelTol',1E-3,...
    'Events',@event_landing); 


settings.m0 = 63.30421;                  % [kg]   Overall Mass
settings.ms = 44.40621;                  % [kg]   Structural Mass (Burnout)
settings.mp = settings.m0-settings.ms;   % [kg]   Propellant Mass
settings.mnc = 6.21;                     % [kg]   Nosecone Mass

% PARACHUTES DETAILS

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

settings.rocket_name = "";
settings.terrain = false;

%% STARTING CONDITIONS

% State
X0 = [0 0 -1000]';
V0 = [0 0 0]';
% V0 = [3;3;0];

%% WIND GENERATION

uw = 0.1:0.1:10;

%% DROGUE 1 + 2
zdrg2 = 500:-20:200;                    % [m] Altitude of drogue 2 opening
Rcirc = 700;
wind_max = zeros(length(zdrg2),1);
R = zeros(length(zdrg2),length(uw));
a = zeros(length(zdrg2),2);
figure

for i = 1:length(zdrg2)
    for j = 1:length(uw)
        para = 1; % Flag for Drogue 1
        X0d1 = [X0;V0];
        settings.zdrg2 = zdrg2(i);
        [Td1,Yd1] = ode113(@descent_easy,[0,1000],X0d1,optionsdrg1,settings,uw(j),0,0,para);
       
        para = 2; % Flag for Drogue 2
        X0d2 = Yd1(end,:);
        
        [Td2,Yd2] = ode113(@descent_easy,[Td1(end),100],X0d2,optionsdrg2,settings,uw(j),0,0,para);
        
        R(i,j) = norm([Yd2(end,1),Yd2(end,2)]);
        
    end
    
    a(i,:) = polyfit(uw,R(i,:),1);
    wind_max(i) = (Rcirc - a(i,2))/a(i,1);
end

plot(zdrg2,wind_max)
xlabel('2nd drogue opening altitude [m] '); ylabel('Max wind [m/s]');



% rogallo wing
% The drogue parachute effects are neglected
settings.para3.S = 15;                           % [m^2]   Surface
settings.para3.mass = 1.466;                     % [kg]   Parachute Mass
settings.para3.CD = 0.4;                         % [/] Parachute Drag Coeff
settings.para3.CL = 0.8;                         % [/] Parachute Lift Coefficient

%% WIND GENERATION

uw = 0.1:0.1:10;
vw = 0;
ww = 0;

%% ROGALLO
zrog = 360:-20:200;                            % [m] Altitude of Rogallo Opening
Rcirc = 700;
wind_max = zeros(length(zrog),1);
R = zeros(length(zrog),length(uw));
a = zeros(length(zrog),2);
h = zeros(length(zrog),1);

figure

for i = 1:length(zrog)
    for j = 1:length(uw)
        X0 = [0 0 -zrog(i)]';
        V0 = [uw(j) 0 0.1]';
        X0d1 = [X0;V0];
        [Td1,Yd1] = ode113(@descent_easy,[0,1000],X0d1,optionsdrg2,settings,uw(j),0,0,3);
        R(i,j) = norm([Yd1(end,1),Yd1(end,2)]);
        
    end
   
    a(i,:) = polyfit(uw,R(i,:),1);
    wind_max(i) = (Rcirc - a(i,2))/a(i,1);
end

plot(zrog,wind_max)
xlabel('rogallo launch altitude [m]'); ylabel('Max wind [m/s]');
        

