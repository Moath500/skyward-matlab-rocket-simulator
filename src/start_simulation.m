<<<<<<< HEAD
% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

close all
clear, clc

%% caricamento dati
run('config.m');

global bool contatore;
global t_plot T_plot alpha_plot beta_plot M_plot CA_plot Drag_plot Forces_plot;
% global WIND alt;


% if bool = 0 the trends during the integration are NOT requested
% if bool = 1 the trends during the integration are saved and plotted
bool = 0;
if bool == 1
    contatore = 1;
    t_plot = 0;
    T_plot = 0;
    alpha_plot = 0;
    beta_plot = 0;
    M_plot = 0;
    CA_plot = 0;
    Drag_plot = 0;
    Forces_plot = 0;
%     WIND = [];
%     alt = [];
end

%% lancio simulatore
[T,Y, Ta,Ya] = MAIN(settings);

if (bool == 1)
    figure();
    plot(t_plot, T_plot, '.'), title('Thrust vs time'), grid on;
    ylabel('Thrust [N]')
    figure();
    plot(t_plot, alpha_plot*180/pi), title('alpha vs time'), grid on;
    ylabel('alpha [deg]')
    figure();
    plot(t_plot, beta_plot*180/pi), title('beta vs time'), grid on;
    ylabel('beta [deg]')

    figure();
    plot(t_plot, M_plot), title('mach vs time'), grid on;
    ylabel('Mach M [-]')

    figure();
    plot(t_plot, CA_plot), title('Aerodyn Coeff vs time'), grid on;
    ylabel('Drag Coeff CD [-]')

    figure();
    plot(t_plot, Drag_plot), title('Drag vs time'), grid on;
    ylabel('Drag D [N]')

    figure();
    plot(t_plot, Forces_plot), title('Axial Force vs time'), grid on;
    ylabel('Axial force [N]')

%     figure()
%     title('Wind'), grid on, hold on;
%     plot(WIND(:, 1), alt)
%     plot(WIND(:, 2), alt)
%     plot(WIND(:, 3), alt)
%     xlabel('Wind magnitude [m/s]')
%     ylabel('Altitude [m]')
%     legend('North','East','Down')
end


[apogee, i_apogee] = max(-Y(:,3));

%% posizione
x = Y(:,1);
y = Y(:,2);
z = -Y(:,3);
X = [x, y, z];

%% velocit�
N = length(x);
%     % derivata centrale
%     vx = (x(3:N)-x(1:N-2))./(T(3:N)-T(1:N-2));
%     vy = (y(3:N)-y(1:N-2))./(T(3:N)-T(1:N-2));
%     vz = (z(3:N)-z(1:N-2))./(T(3:N)-T(1:N-2));
%     % aggiungo derivata estremi
%     vx = [0; vx; (x(end)-x(end-1))/(T(end)-T(end-1))];
%     vy = [0; vy; (y(end)-y(end-1))/(T(end)-T(end-1))];
%     vz = [0; vz; (z(end)-z(end-1))/(T(end)-T(end-1))];
vx = Y(:,4);
vy = Y(:,5);
vz = -Y(:,6);
V = [vx, vy, vz];

%% accelerazione
% derivata centrale
ax = (vx(3:N)-vx(1:N-2))./(T(3:N)-T(1:N-2));
ay = (vy(3:N)-vy(1:N-2))./(T(3:N)-T(1:N-2));
az = (vz(3:N)-vz(1:N-2))./(T(3:N)-T(1:N-2));
% aggiungo derivata estremi
ax = [vx(2)/T(2); ax; (vx(end)-vx(end-1))/(T(end)-T(end-1))];
ay = [vy(2)/T(2); ay; (vy(end)-vy(end-1))/(T(end)-T(end-1))];
az = [vz(2)/T(2); az; (vz(end)-vz(end-1))/(T(end)-T(end-1))];

% % derivata in avanti
% ax = (vx(2:N)-vx(1:N-1))./(T(2:N)-T(1:N-1));
% ay = (vy(2:N)-vy(1:N-1))./(T(2:N)-T(1:N-1));
% az = (vz(2:N)-vz(1:N-1))./(T(2:N)-T(1:N-1));
% % aggiungo derivata estremi
% ax = [0; ax];
% ay = [0; ay];
% az = [0; az];

A = [ax, ay, az];

clear ax ay az
clear vx vy

mod_X = zeros(length(X), 1);
mod_V = mod_X;
mod_A = mod_X;

for k = 1:length(X)
    mod_X(k) = norm(X(k,:));
    mod_V(k) = norm(V(k,:));
    mod_A(k) = norm(A(k,:));
end
[max_dist, imax_dist] = max(mod_X);
[max_v, imax_v] = max(mod_V);
[max_a, imax_a] = max(mod_A);

M = zeros(length(X), 1);
Tamb = M;
Ttot = Tamb;
for n = 1:length(M)
    h = z(n);
    if (h < 0)
        h = 0;
    end
    [Tamb(n), a, ~, ~] = atmoscoesa(h);
    M(n) = mod_V(n)/a;
    Ttot(n) = Tamb(n)*(1+0.2*M(n)^2);
end

[max_M, imax_M] = max(M);

%% launch pad exit index
iexit = find(mod_X <= settings.lrampa);
iexit = iexit(end);


%% DATA RECORD (display)
disp(' ')
disp('DATA RECORD:')
fprintf('apogee reached: %g [m] \n', apogee);
fprintf('apogee velocity: %g [m/s] \n', mod_V(i_apogee));
fprintf('time: %g [sec] \n\n', Ta(end))

fprintf('max speed reached: %g [m/s] \n', max_v)
fprintf('altitude: %g [m] \n', z(imax_v))
fprintf('Mach: %g [-] \n', M(imax_v))
fprintf('time: %g [sec] \n\n', T(imax_v))

fprintf('max Mach reached: %g [-] \n', max_M)
fprintf('altitude: %g [m] \n', z(imax_M))
fprintf('veocity: %g [m/s] \n', mod_V(imax_M))
fprintf('time: %g [sec] \n\n', T(imax_M))

fprintf('max acceleration reached: %g [m/s2] = %g [g] \n', max_a, max_a/9.80665)
fprintf('veocity: %g [m/s] \n', mod_V(imax_a))
fprintf('time: %g [sec] \n\n', T(imax_a))

fprintf('run on launch pad: %g [m] \n', mod_X(iexit))
fprintf('speed at launch pad exit: %g [m/s] \n', mod_V(iexit))
fprintf('time: %g [sec] \n\n', T(iexit))

%% plot some interesting stuff
plots = true;
if plots
    % figure(), hold on, grid on;
    % plot(T, z);
    % plot(T, vz);
    % plot(T, mod_V);
    % plot(T, mod_A/9.80665)
    % legend('z', 'vz', 'mod V', 'mod A');
    % xlabel('Time [sec]');
    % ylabel('z [m], vz [m/s], mod V [m/s], mod A [g]')

    figure();
    subplot(3,1,1)
    plot(T, z), grid on, xlabel('time [sec]'), ylabel('altitude [m]');
    subplot(3,1,2)
    plot(T, vz), grid on, hold on;
    plot(T, mod_V);
    xlabel('time [sec]'), ylabel('vz [m/s], |V| [m/s]');
    legend('vz', '|V|');
    subplot(3,1,3)
    plot(T, mod_A/9.80665), grid on;
    xlabel('time [sec]'), ylabel('|A| [g]');
    legend('|A|');

    figure();
    plot(T, z), grid on, xlabel('time [sec]'), ylabel('altitude [m]');

    figure();
    plot(T, vz), grid on, hold on;
    plot(T, mod_V);
    xlabel('time [sec]'), ylabel('vz [m/s], |V| [m/s]');
    legend('vz', '|V|');

    figure();
    plot(T, mod_A/9.80665), grid on;
    xlabel('time [sec]'), ylabel('|A| [g]');
    legend('|A|');

    figure();
    plot3(y, x, z), axis equal, hold on, grid on;
    plot3(Ya(end,2), Ya(end,1), -Ya(end,3), '*')
    plot3(0, 0, 0, 'or')
    % visualizzazione circonferenze distanze
    theta_plot = linspace(0,2*pi);
    R_plot = [1, 2, 3, 4, 5]*1000;
    for j = 1:length(R_plot)
        x_plot = R_plot(j)*cos(theta_plot');
        y_plot = R_plot(j)*sin(theta_plot');
        z_plot = zeros(length(theta_plot), 1);
        plot3(y_plot, x_plot, z_plot, '--r')
    end

    title('Trajectory')
    xlabel('y, East [m]'), ylabel('x, North [m]'), zlabel('Altitude [m]')

    figure();
    plot(y, x), axis equal, hold on, grid on;
    % visualizzazione circonferenze distanze
    theta_plot = linspace(0,2*pi);
    R_plot = [1, 2, 3, 4, 5]*1000;
    for j = 1:length(R_plot)
        x_plot = R_plot(j)*cos(theta_plot');
        y_plot = R_plot(j)*sin(theta_plot');
        plot(y_plot, x_plot, '--r')
    end
    xlabel('y, East [m]'), ylabel('x, North [m]');

    % lengths to 0 if they are less than 1 meter
    x_flight=x;
    y_flight=y;
    z_flight=z;
    for ii = 1 : length(x)
        if norm(x_flight(ii))<1
            x_flight(ii) = 0;
        end
        if norm(y_flight(ii))<1
            y_flight(ii) = 0;
        end
        if norm(z_flight(ii))<1
            z_flight(ii) = 0;
        end
    end

    figure();
    subplot(1,3,1);
    plot(y_flight/1000, x_flight/1000), axis equal, hold on, grid on;
    plot(Ya(end,2)/1000, Ya(end,1)/1000, '*');
    plot(0, 0, 'or')
    xlabel('y, East [Km]'), ylabel('x, North [Km]');
    subplot(1,3,2);
    plot(x_flight/1000, z_flight/1000), hold on, grid on;
    plot(Ya(end,1)/1000, -Ya(end,3)/1000, '*');
    plot(0, 0, 'or')
    xlabel('x, North [Km]'), ylabel('z, Altitude [Km]');
    %setting limit if is parallel to east
    if sum(x_flight) == 0
        xlim([-1 1]);
    end
    subplot(1,3,3);
    plot(y_flight/1000, z_flight/1000), hold on, grid on;
    plot(Ya(end,2)/1000, -Ya(end,3)/1000, '*');
    plot(0, 0, 'or')
    xlabel('y, East [Km]'), ylabel('z, Altitude [Km]');
    %setting limit if is parallel to north
    if sum(y_flight) == 0
        xlim([-1 1]);
    end

    % angular rates
    figure();
    plot(Ta, Ya(:, 7)*180/pi)
    hold on, grid on
    plot(Ta, Ya(:, 8)*180/pi)
    plot(Ta, Ya(:, 9)*180/pi)
    xlabel('time [sec]'), ylabel('p, q, r [grad/sec]')
    legend('p: roll rate', 'q: pitch rate', 'r: yaw rate')
    title('angular rates vs time')

    figure();
    subplot(3,1,1)
    plot(Ta, Ya(:, 8)*180/pi), grid on;
    xlabel('time [sec]'), ylabel('pitch rate q [grad/sec]')
    subplot(3,1,2)
    plot(Ta, Ya(:, 9)*180/pi), grid on;
    xlabel('time [sec]'), ylabel('yaw rate r [grad/sec]')
    subplot(3,1,3)
    plot(Ta, Ya(:, 7)*180/pi), grid on;
    xlabel('time [sec]'), ylabel('roll rate p [grad/sec]')


    % angle
    alpha = zeros(length(Ta),1);
    beta = zeros(length(Ta),1);
    roll = zeros(length(Ta),1);

    for k = 2:length(Ta)
        alpha(k) = alpha(k-1) + (Ya(k, 8) + Ya(k-1, 8))/2*180/pi*(T(k)-T(k-1));
        beta(k) = beta(k-1) + (Ya(k, 9) + Ya(k-1, 9))/2*180/pi*(T(k)-T(k-1));
        roll(k) = roll(k-1) + (Ya(k, 7) + Ya(k-1, 7))/2*180/pi*(T(k)-T(k-1));
    end
    figure(), plot(Ta, alpha+settings.OMEGA*180/pi), title('Alpha vs time'), grid on;
    figure(), plot(Ta, beta), title('Beta vs time'), grid on;
    figure(), plot(Ta, roll), title('Roll vs time'), grid on;


    % Mach
    figure(), plot(T, M), grid on;
    xlabel('time [sec]'), ylabel('Mach [-]');

    % Stagnation Temperature
    figure();
    plot(T, Tamb-273.15);
    hold on, grid on;
    plot(T, Ttot-273.15)
    title('Temperature profile')
    xlabel('time [sec]'), ylabel('Temperature [�C]');
    legend('Surrounding', 'Total', 'location', 'best')

end
