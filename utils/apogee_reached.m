<<<<<<< HEAD
% TARGET ??
% usare il simulatore a 6dof per simulare diverse configuarzioni del
% missile e vedere la quota massima raggiunta per poter fare un confronto e
% decidere la soluzione pi� perormante

%% COSA MI SERVE ??
% 1) i vari file di DATCOM per le diverse configurazioni (geometrie e CG)
% 2) le inerzie del missile per le diverse configurazioni
% 3) il profilo della spinta (rocket engine data)

%% COME OTTENERE LA QUOTA DELL'APOGEO DAL SIMULATORE ??

% ODE-Function of the 6DOF Rigid Rocket Model
% State = ( x y z | u v w | p q r | q0 q1 q2 q3 )
%
% (x y z): NED Earth's Surface Centered Frame ("Inertial") coordinates
% (u v w): body frame velocities
% (p q r): body frame angular rates
% (q0 q1 q2 q3): attitude unit quaternion

% % [T,Y]=MAIN();
% %
% % apogeo = max(-Y(:,3));
% % fprintf('Apogee: %g [m] \n', apogeo);

%% CODICE
% close all
clear, clc

global bool
bool = 0;
% if bool = 0 the trends during the integration are NOT requested
% if bool = 1 the trends during the integration are saved and plotted


% se la geometria non cambia faccio caricare i dati con i coeff
% aerodinamici da "config.m" solo una volta e poi cambio le inerzie dentro
% al ciclo "for"
run('config.m');

% masse
m_fuel = 18.6; %kg
% m_empty = 40.4;
% m_empty = [44.6:]; %kg
% m_full = m_empty + m_fuel;
m_full = linspace(44.6, 84.6, 30)';
m_empty = m_full - m_fuel;


% dimensioni missile
D = 0.174; %m
% D = 0.214; %m
R = D/2; %m
S = pi*D^2/4; %m2
L = 4.30; %m

% inerzie
Ixx_full = m_full.*R^2/2;
Ixx_empty = m_empty.*R^2/2;
Iyy_full = m_full.*(R^2/4 + L^2/3);
Iyy_empty = m_empty.*(R^2/4 + L^2/3);
Izz_full = Iyy_full;
Izz_empty = Iyy_empty;

n_casi = length(Ixx_full);

apogeo = zeros(n_casi, 1);
max_dist = zeros(n_casi, 1);
max_v = zeros(n_casi, 1);
max_a = zeros(n_casi, 1);

for n = 1:n_casi
    % carica i dati dai file di DATCOM per ogni caso e salvali nel file
    % .mat per il run del simulatore

%     % se invece per ogni caso cambia anche la geometria
%     file_name = ['for006_empty(',num2str(n),').mat'];
%     load(file_name);
%     save('for006_empty.mat', 'Coeffs', 'State');
%
%     file_name = ['for006_full(',num2str(n),').mat'];
%     save('for006_full.mat', 'Coeffs', 'State');
%
%     run('config.m');

    % setto la massa e le inerzie del caso n-th
    settings.ms = m_empty(n);
    settings.m0 = m_full(n);

    settings.Ixxf = Ixx_full(n); %Inertia to x-axis (Full)
    settings.Ixxe = Ixx_empty(n); %Inertia to x-axis (Empty)
    settings.Iyyf = Iyy_full(n); %Inertia to y-axis (Full)
    settings.Iyye = Iyy_empty(n); %Inertia to y-axis (Empty)
    settings.Izzf = Izz_full(n); %Inertia to z-axis (Full)
    settings.Izze = Izz_empty(n); %Inertia to z-axis (Empty)

%     % prova senza inerzie
%     settings.Ixxf=2.1; %Inertia to x-axis (Full)
%     settings.Ixxe=1.4; %Inertia to x-axis (Empty)
%     settings.Iyyf=20; %Inertia to y-axis (Full)
%     settings.Iyye=15; %Inertia to y-axis (Empty)
%     settings.Izzf=20; %Inertia to z-axis (Full)
%     settings.Izze=15; %Inertia to z-axis (Empty)

%     save('settings.mat', 'settings')

    % faccio girare il simulatore
    [T,Y, Ta,Ya] = MAIN(settings);

    [apogeo(n), iapogeo] = max(-Y(:,3));
    fprintf('apogee is equal to %g m \n', apogeo(n));
    fprintf('apogeo raggiunto in %g sec \n\n', Ta(end))

%     % plot vari
%     % posizione
%     x = Y(:,1);
%     y = Y(:,2);
%     z = -Y(:,3);
%     X = [x, y, z];
%
%     % velocit�
%     N = length(x);
% %     % derivata centrale
% %     vx = (x(3:N)-x(1:N-2))./(T(3:N)-T(1:N-2));
% %     vy = (y(3:N)-y(1:N-2))./(T(3:N)-T(1:N-2));
% %     vz = (z(3:N)-z(1:N-2))./(T(3:N)-T(1:N-2));
% %     % aggiungo derivata estremi
% %     vx = [0; vx; (x(end)-x(end-1))/(T(end)-T(end-1))];
% %     vy = [0; vy; (y(end)-y(end-1))/(T(end)-T(end-1))];
% %     vz = [0; vz; (z(end)-z(end-1))/(T(end)-T(end-1))];
%     vx = Y(:,4);
%     vy = Y(:,5);
%     vz = -Y(:,6);
%     V = [vx, vy, vz];
%
%     % accelerazione
%     % derivata centrale
%     ax = (vx(3:N)-vx(1:N-2))./(T(3:N)-T(1:N-2));
%     ay = (vy(3:N)-vy(1:N-2))./(T(3:N)-T(1:N-2));
%     az = (vz(3:N)-vz(1:N-2))./(T(3:N)-T(1:N-2));
%     % aggiungo derivata estremi
%     ax = [vx(2)/T(2); ax; (vx(end)-vx(end-1))/(T(end)-T(end-1))];
%     ay = [vy(2)/T(2); ay; (vy(end)-vy(end-1))/(T(end)-T(end-1))];
%     az = [vz(2)/T(2); az; (vz(end)-vz(end-1))/(T(end)-T(end-1))];
%     A = [ax, ay, az];
%
%     mod_X = zeros(length(X), 1);
%     mod_V = mod_X;
%     mod_A = mod_X;
%
%     for k = 1:length(X)
%         mod_X(k) = norm(X(k,:));
%         mod_V(k) = norm(V(k,:));
%         mod_A(k) = norm(A(k,:));
%     end
%     [max_dist, imax_dist] = max(mod_X);
%     [max_v, imax_v] = max(mod_V);
%
%     % plot
%     figure(), hold on, grid on;
%     plot(T, z);
%     plot(T, vz);
%     plot(T, mod_V);
%     legend('z', 'vz', 'mod V');
end

% istogramma apogeo vs massa struttura
figure(), grid on;
stem(m_full, apogeo);
xlabel('overall mass = structural mass + fuel(=18.6kg) [kg]');
ylabel('apogee [m]')

[maxApogeo, imaxApogeo] = max(apogeo);
fprintf('massimo apogeo raggiunto: %g m \n', maxApogeo);
fprintf('ottenuto con massa strutturale di %g kg \n\n', m_empty(imaxApogeo));
