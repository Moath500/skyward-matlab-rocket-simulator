% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

function [apogee, max_a,Vexit,t,vect_XCP]=start_simulation(settings)


path = genpath(pwd);
addpath(path);


%% START THE CHOSEN SIMULATION
% T = vector of time used by ODE, [s] also for Tf Ta
% Y = State = ( x y z | u v w | p q r | q0 q1 q2 q3 ) also for Ya,Yf corresponding to T

        [T,Y,Ta,Ya,bound_value] = std_run_ballistic(settings);

%% NOT-STOCHASTIC SIMULATIONS (N=1)

    N = length(Y(:,1));
    n = length(Ya(:,1));
    
    
    % POSITIONS
    
    x = Y(:,1);
    y = Y(:,2);
    z = -Y(:,3);
    X = [x, y, z];
    
    z_a  = z(1:n);
    
    [apogee, i_apogee] = max(-Y(:,3)); % position, index of position at apogee
    
    
    % VELOCITIES
    
    
    u = Y(:,4);
    v = Y(:,5);
    w = -Y(:,6);
    V = [u, v, w];
    
    % ACCELERATIONS
    
    % main derivatives
    ax = (u(3:N)-u(1:N-2))./(T(3:N)-T(1:N-2));
    ay = (v(3:N)-v(1:N-2))./(T(3:N)-T(1:N-2));
    az = (w(3:N)-w(1:N-2))./(T(3:N)-T(1:N-2));
    
    % add derivative at the boundaries
    ax = [u(2)/T(2); ax; (u(end)-u(end-1))/(T(end)-T(end-1))];
    ay = [v(2)/T(2); ay; (v(end)-v(end-1))/(T(end)-T(end-1))];
    az = [w(2)/T(2); az; (w(end)-w(end-1))/(T(end)-T(end-1))];
    
    A = [ax, ay, az];
    
    clear('ax', 'ay', 'az', 'vx', 'vy')
    
    % MAXIMUM POSITIONS, VELOCITIES AND ACCELERATIONS
    
    % pre-allocation
    abs_X = zeros(length(X), 1);
    abs_V = abs_X;
    abs_A = abs_X;
    
    % determine the module of every row element
    for k = 1:length(X)
        abs_X(k) = norm(X(k,:));
        abs_V(k) = norm(V(k,:));
        abs_A(k) = norm(A(k,:));
    end
    
    abs_Va = abs_V(1:n);
    abs_Aa = abs_A(1:n);
    
    
    [max_dist, imax_dist] = max(abs_X);
    [max_v, imax_v] = max(abs_V);
    [max_a, imax_a] = max(abs_A);
    iexit = find(abs_X <= settings.lrampa);  % checking where the missile is undocked from the hook of the launch pad
    iexit = iexit(end);
    
    % TEMPERATURE AND MACH NUMBER
    
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
        M(n) = abs_V(n)/a;
        Ttot(n) = Tamb(n)*(1+0.2*M(n)^2);
    end
    % determine the maximum Mach number
    [max_M, imax_M] = max(M);
    
    load data_ascent.mat
    delete data_ascent.mat
    
    
    max_a=max_a/9.80665;
    Vexit=abs_V(iexit);
    vect_XCP=-data_ascent.coeff.XCP;
    t=data_ascent.integration.t;
    
