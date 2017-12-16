% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

close all
clear all
clc

%% LOAD DATA

run('config.m');

%% START THE CHOSEN SIMULATION
% T = vector of time used by ODE, [s] also for Tf Ta
% Y = State = ( x y z | u v w | p q r | q0 q1 q2 q3 ) also for Ya,Yf corresponding to T  

tic

% Checking if stochastic or standard simulation needed
if settings.ballistic
    if settings.stoch.N > 1
        fprintf('Stochastic Ballistic Simulation Started...\n\n');
        if settings.stoch.parallel
            [LP,Z] = stoch_run_bal_p(settings);
        else
            [LP,Z] = stoch_run_bal(settings);
        end
    else
        fprintf('Standard Ballistic Simulation Started...\n\n');
        [T,Y, Ta,Ya] = std_run_ballistic(settings);
    end
else
    if settings.stoch.N > 1
        fprintf('Stochastic Simulation Started...\n\n');
        if settings.stoch.parallel
            [LP,Z] = stoch_run_p(settings);
        else
            [LP,Z] = stoch_run(settings);
        end
    else
        fprintf('Standard Simulation Started...\n\n');
        [T,Y, Ta,Ya] = std_run(settings);
    end
end
toc

if settings.stoch.N == 1

%% POSITIONS

x = Y(:,1);
y = Y(:,2);
z = -Y(:,3);
X = [x, y, z];
[apogee, i_apogee] = max(-Y(:,3)); % position, index of position at apogee

%% VELOCITIES

N = length(x);
vx = Y(:,4);
vy = Y(:,5);
vz = -Y(:,6);
V = [vx, vy, vz];

%% ACCELERATIONS

% main derivatives
ax = (vx(3:N)-vx(1:N-2))./(T(3:N)-T(1:N-2));
ay = (vy(3:N)-vy(1:N-2))./(T(3:N)-T(1:N-2));
az = (vz(3:N)-vz(1:N-2))./(T(3:N)-T(1:N-2));

% add derivative at the boundaries
ax = [vx(2)/T(2); ax; (vx(end)-vx(end-1))/(T(end)-T(end-1))];
ay = [vy(2)/T(2); ay; (vy(end)-vy(end-1))/(T(end)-T(end-1))];
az = [vz(2)/T(2); az; (vz(end)-vz(end-1))/(T(end)-T(end-1))];

A = [ax, ay, az];

clear('ax', 'ay', 'az', 'vx', 'vy')

%% MAXIMUM POSITIONS, VELOCITIES AND ACCELERATIONS

% pre-allocation 
mod_X = zeros(length(X), 1);
mod_V = mod_X;
mod_A = mod_X;

% determine the module of every row element
for k = 1:length(X)
    mod_X(k) = norm(X(k,:));
    mod_V(k) = norm(V(k,:));
    mod_A(k) = norm(A(k,:));
end

[max_dist, imax_dist] = max(mod_X);
[max_v, imax_v] = max(mod_V);
[max_a, imax_a] = max(mod_A);
iexit = find(mod_X <= settings.lrampa);  % checking where the missile is undocked from the hook of the launch pad
iexit = iexit(end);

%% TEMPERATURE AND MACH NUMBER

% pre-allocation
M = zeros(length(X), 1);
Tamb = M;
Ttot = Tamb;
for n = 1:length(M)
    h = z(n);
    if (h < 0)
        h = 0;
    end
    [Tamb(n), a, ~, ~] = atmosisa(h);
    M(n) = mod_V(n)/a;
    Ttot(n) = Tamb(n)*(1+0.2*M(n)^2);
end
% determine the maximum Mach number
[max_M, imax_M] = max(M);

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

%% PLOT

if settings.standard_plot
    run('standard_plot.m')
end

if settings.ascend_plot
    run('ascend_plot.m') 
end

clear('ascend_plot.mat')

end
