function [ dY ] = descent_parachute( t,Y,settings,uw,vw,ww,para )
% ODEFUN for Parachute descent
% State = ( x y z | u v w  )

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

x = Y(1);
y = Y(2);
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);


% global bool
% % if bool = 0 the trends during the integration are NOT requested
% % if bool = 1 the trends during the integration are saved and plotted


% Adding Wind (supposed to be added in NED axes);

% constant wind
wind = [uw vw ww];

% % Wind model
% h = 0;
% if -z > 0
%     h = -z+settings.z0;
% end
% 
% % global alt
% % if bool == 1
% %     alt = [alt; h];
% % end
% 
% % wind in NED axis
% [uw,vw] = atmoshwm07(settings.wind.Lat, settings.wind.Long, h, ...
%     'day',settings.wind.Day,'model','quiet');
% wind = [uw, vw, ww];
% 
% % global WIND
% % if bool == 1
% %     WIND = [WIND; uw, vw, ww];
% % end

% % Wind Model
% [uw,vw,ww] = wind_generator(settings,z,t,Q);
% wind = [ uw,vw,ww ];

% Adding wind;
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

V_norm = norm([ur vr wr]);

% Constants
switch para
    case 1
        S = settings.para1.S;   %Parachute Surface
        CD = settings.para1.CD; %Parachute CD
        CL = settings.para1.CL; %Parachute CL
        pmass = 0;              %detached parachute mass
    case 2
        S = settings.para2.S;
        CD = settings.para2.CD;
        CL = settings.para2.CL;
        pmass = settings.para1.mass;
    case 3
        S = settings.para3.S;
        CD = settings.para3.CD;
        CL = settings.para3.CL;
        pmass = settings.para1.mass + settings.para2.mass;
    otherwise
end

rad = sqrt(S/pi); %Equivalent parachute radius;

%Side Surface 
Ss = pi*(rad/2)^2;

g = 9.81;

m = settings.ms - pmass;

%Atmosphere
if -z<0
    z = 0;
end

[~, ~, ~, rho] = atmoscoesa(-z);

%Center of Mass



%Body Frame
D = 0.5*rho*V_norm^2*S*CD;
L = 0.5*rho*V_norm^2*S*CL;

% HP: il paracadute è approssimato come una superficie rettangolare
% orientata in modo da essere sempre normale alla direzione della velocita
% relativa all'aria e sempre con un asse orizzontale e la forza Lift verso
% l'alto
% in modo da puntare ad atterrare nel punto [x,y] = (0m a nord,500m ad
% ovest) rispetto la postazione di lancio


tt = [ur vr wr]; % versore tangenziale
kk = [-vr ur 0]; % versore orizzontale
if (kk(1) == 0) && (kk(2) == 0)
    kk = [-v u 0];
end
tt = tt/norm(tt);
kk = -kk/norm(kk);
nn = cross(tt, kk); nn = nn/norm(nn); % versore normale
if (nn(3) > 0) % nn diretto verso il basso ?
    nn = cross(kk, tt); nn = nn/norm(nn);
end


%Forces

% X = D*ur/V_norm;
% Y = D*vr/V_norm;
% Z = D*wr/V_norm -m*g;

% F = [-X,-Y,-Z]';
F = -D*tt' + L*nn' + m*g*[0 0 1]';


du=F(1)/m;
dv=F(2)/m;
dw=F(3)/m;





dY(1:3) = [u v w]';
dY(4) = du;
dY(5) = dv;
dY(6) = dw;

dY = dY';

end