% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016


    %% ALTITUDE
    
    figure();
    plot(T, z), grid on, xlabel('time [s]'), ylabel('altitude [m]'), title('Altitude vs time');
    
    %% VELOCITY
    
    figure();
    plot(T, vz), grid on, hold on;
    plot(T, mod_V);
    xlabel('time [s]'), ylabel('vz [m/s], |V| [m/s]'), title('Velocities vs time');
    legend('vz', '|V|');
    
    %% ACCELERATION
    
    figure();
    plot(T, mod_A/9.80665), grid on;
    xlabel('time [s]'), ylabel('|A| [g]'), title('Acceleration vs time');
    legend('|A|');
    
    %% ALTITUDE,VELOCITY,ACCELERATION(subplotted)
    
    figure()
    subplot(3,1,1)
    plot(T, z), grid on, xlabel('time [s]'), ylabel('altitude [m]'), title('Altitude vs time');
    
    subplot(3,1,2)
    plot(T, vz), grid on, hold on;
    plot(T, mod_V);
    xlabel('time [s]'), ylabel('vz [m/s], |V| [m/s]');
    legend('vz', '|V|'), title('Velocities vs time');
    
    subplot(3,1,3)
    plot(T, mod_A/9.80665), grid on;
    xlabel('time [s]'), ylabel('|A| [g]');
    legend('|A|');
    title('Acceleration vs time')
    
    %% 3D TRAJECTORY
    
    figure();
    plot3(y, x, z), axis equal, hold on, grid on;
    plot3(Ya(end,2), Ya(end,1), -Ya(end,3), '*')
    plot3(0, 0, 0, 'or')
    % concentric circles distance 
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
    
    %% TRAECTORY'S TOP VIEW
    
    figure();
    plot(y, x), axis equal, hold on, grid on;
    p1 = plot(0,0,'r.','MarkerSize',21);
    p2 = plot(y(end),x(end),'rx','MarkerSize',14);
    title('Trajectory Top view')
    xlabel('y, East [m]'), ylabel('x, North [m]');
    legend([p1,p2],{'Launch Point','Landing Point'})
    clear p1 p2
    
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
    
    %% TRAJECTORY PROJECTIONS(subplotted)
    
    figure()
    suptitle('Trajectory projections')
    subplot(1,3,1)
    plot(y_flight/1000, x_flight/1000), axis equal, hold on, grid on;
    p11 = plot(Ya(end,2)/1000, Ya(end,1)/1000, 'r*','markersize',7);
    p12 = plot(0, 0, 'r.','markersize',14);
    p13 = plot(Y(end,2)/1000, Y(end,1)/1000, 'rx','markersize',7);
    xlabel('y, East [Km]'), ylabel('x, North [Km]');
    
    subplot(1,3,2);
    plot(x_flight/1000, z_flight/1000), hold on, grid on;
    plot(Ya(end,1)/1000, -Ya(end,3)/1000, 'r*','markersize',7);
    plot(Y(end,1)/1000, -Y(end,3)/1000, 'rx','markersize',7);
    plot(0, 0, 'r.','markersize',14);
    xlabel('x, North [Km]'), ylabel('z, Altitude [Km]');
    %setting limit if is parallel to east
    if sum(x_flight) == 0
        xlim([-1 1]);
    end
    
    subplot(1,3,3);
    plot(y_flight/1000, z_flight/1000), hold on, grid on;
    plot(Ya(end,2)/1000, -Ya(end,3)/1000, 'r*','markersize',7);
    plot(Y(end,2)/1000, -Y(end,3)/1000, 'rx','markersize',7);
    plot(0, 0, 'r.','markersize',14);
    xlabel('y, East [Km]'), ylabel('z, Altitude [Km]');
    %setting limit if is parallel to north
    if sum(y_flight) == 0
        xlim([-1 1]);
    end
    legend([p11, p12, p13],{'Apogee','Launch Point','Landing Point'})
    
    %% ANGULAR RATES (combined)
    
    figure();
    plot(Ta, Ya(:, 7)*180/pi)
    hold on, grid on
    plot(Ta, Ya(:, 8)*180/pi)
    plot(Ta, Ya(:, 9)*180/pi)
    xlabel('time [s]'), ylabel('p, q, r [grad/s]')
    legend('p: roll rate', 'q: pitch rate', 'r: yaw rate')
    title('Angular rates vs time')
    
    %% ANGULAR RATES (subplotted)
    
    figure(),suptitle('Angular rates vs time');
    subplot(3,1,1)
    plot(Ta, Ya(:, 8)*180/pi), grid on;
    xlabel('time [s]'), ylabel('pitch rate q [grad/s]')
    
    subplot(3,1,2)
    plot(Ta, Ya(:, 9)*180/pi), grid on;
    xlabel('time [s]'), ylabel('yaw rate r [grad/s]')
    
    subplot(3,1,3)
    plot(Ta, Ya(:, 7)*180/pi), grid on;
    xlabel('time [s]'), ylabel('roll rate p [grad/s]')
    
    
    %% EULERIAN ANGLES (3 figures)
    
    pitch_angle = zeros(length(Ta),1);
    yaw_angle = zeros(length(Ta),1);
    roll_angle = zeros(length(Ta),1);
    
    for k = 2:length(Ta)
        pitch_angle(k) = pitch_angle(k-1) + (Ya(k, 8) + Ya(k-1, 8))/2*180/pi*(T(k)-T(k-1));
        yaw_angle(k) = yaw_angle(k-1) + (Ya(k, 9) + Ya(k-1, 9))/2*180/pi*(T(k)-T(k-1));
        roll_angle(k) = roll_angle(k-1) + (Ya(k, 7) + Ya(k-1, 7))/2*180/pi*(T(k)-T(k-1));
    end
    
    figure() 
    plot(Ta, pitch_angle+settings.OMEGA*180/pi)
    title('Pitch angle vs time')
    grid on, xlabel('time [s]'), ylabel('pitch angle [deg]');
    
    figure()
    plot(Ta, yaw_angle)
    title('Yaw angle vs time')
    grid on, xlabel('time [s]'), ylabel('yaw angle [deg]');
    
    figure()
    plot(Ta, roll_angle)
    title('Roll angle vs time')
    grid on, xlabel('time [s]'), ylabel('roll angle [deg]')
    
    
    %% MACH NUMBER
    
    figure(), plot(T, M), grid on;
    xlabel('time [s]'), ylabel('Mach [/]'), title('Mach number vs time');
    
    %% STAGNATION TEMPERATURE
    
    figure();
    plot(T, Tamb-273.15);
    hold on, grid on;
    plot(T, Ttot-273.15)
    title('Temperature profile')
    xlabel('time [s]'), ylabel('Temperature [C]');
    legend('Surrounding', 'Total', 'location', 'best')
    
    clear variable
    