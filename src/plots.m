%{
 Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 13/01/2018

%}

if settings.stoch.N == 1
    
    load('ascent_plot.mat')
    if settings.ballistic
        load('descent_plot.mat')
        delete('descent_plot.mat')
    else
        load('descent_para_plot.mat')
        delete('descent_para_plot.mat')
    end
    
    %% ASCENT PLOTS
    
    figure('Name','Stability Margin - ascent Phase','NumberTitle','off');
    plot(data_ascent.integration.t, -data_ascent.coeff.XCP,'.'), title('Stability margin vs time'), grid on;
    xlabel('Time [s]'); ylabel('S.M.[/]')
    
    figure('Name','Forces - ascent Phase','NumberTitle','off');
    subplot(2,2,1)
    plot(data_ascent.integration.t, data_ascent.forces.T, '.'), grid on;
    xlabel('Time [s]'); ylabel('Thrust [N]');
    
    subplot(2,2,2)
    plot(data_ascent.integration.t, data_ascent.forces.AeroDyn_Forces(:,1)),grid on;
    xlabel('Time [s]'); ylabel('X-body force [N]')
    
    subplot(2,2,3)
    plot(data_ascent.integration.t, data_ascent.forces.AeroDyn_Forces(:,2)), grid on;
    xlabel('Time [s]'); ylabel('Y-body force [N]')
    
    subplot(2,2,4)
    plot(data_ascent.integration.t, data_ascent.forces.AeroDyn_Forces(:,3)), grid on;
    xlabel('Time [s]'); ylabel('Z-body force [N]')
    
    figure('Name','Aerodynamics Angles - ascent Ahase','NumberTitle','off');
    subplot(2,1,1)
    plot(data_ascent.integration.t, data_ascent.interp.alpha*180/pi), grid on;
    xlabel('Time [s]'); ylabel('alpha [deg]')
    
    subplot(2,1,2)
    plot(data_ascent.integration.t, data_ascent.interp.beta*180/pi), grid on;
    xlabel('Time [s]'); ylabel('beta [deg]')
    
    figure('Name','Drag Coefficient - ascent Phase','NumberTitle','off');
    plot(data_ascent.integration.t, data_ascent.coeff.CA), title('Drag Coefficient vs time'), grid on;
    xlabel('Time [s]'); ylabel('Drag Coeff CD [/]')
    
    
    %% 3D TRAJECTORY
    figure('Name','3D Trajectory - All Flight','NumberTitle','off');
    plot3(y,x,z,'Linewidth',2.5), axis equal, hold on, grid on;
    title('Trajectory')
    xlabel('y, East [m]'), ylabel('x, North [m]'), zlabel('Altitude [m]')
    
    % randomly generation of colors:
    Np = settings.Npara;
    Colors = rand(3, Np);
    
    % adding concentric circles
    if not(settings.terrain)
        theta_plot = linspace(0,2*pi);
        R_plot = [1, 2, 3, 4, 5]*1000;
        
        for j = 1:length(R_plot)
            x_plot = R_plot(j)*cos(theta_plot');
            y_plot = R_plot(j)*sin(theta_plot');
            z_plot = zeros(length(theta_plot), 1);
            plot3(y_plot, x_plot, z_plot, '--r')
        end
        
    else
        % adding the terrain map
        X_t = -6000:30:6000;
        Y_t = -6000:30:6000;
        L_X = length(X_t);
        L_Y = length(Y_t);
        Z_t = zeros(L_X,L_Y);
        for i = 1:L_X
            for j = 1:L_Y
                Z_t(i,j) = settings.funZ(X_t(i),Y_t(j));
            end
        end
        surf(Y_t,X_t,-Z_t,'EdgeColor','none'); colorbar; view(0,25);
    end
    
    % adding boundary values
    [~, i_tb] = min(abs(settings.tb-Ta));
    
    h = zeros(Np+3, 1);
    for i = 1:Np
        h(i) = plot3(bound_value(i).X(1), bound_value(i).X(2), bound_value(i).X(3), 'o',...
            'MarkerSize', 10, 'MarkerFaceColor', Colors(:,i) , 'MarkerEdgeColor', 'none');
    end
    
    h(Np+1) = plot3(Yf(end,2),Yf(end,1),-Yf(end,3),'rx','markersize',10);
    h(Np+2) = plot3(0, 0, 0, '*');
    h(Np+3) = plot3(Yf(i_tb,2),Yf(i_tb,1),-Yf(i_tb,3),'ro','MarkerSize',10);
    
    if not(settings.ballistic)
        legend(h(:), {'Apogee', strcat('parachute ',  " " , string(2:Np), " ", 'opening'), 'Landing point',...
            'Launch point', 'burning time position'}, 'Location', 'northeast', 'Fontsize', 18);
    else
        legend(h(:), {'Apogee', 'Landing point', 'Launch point', 'burning time position'},...
            'Location', 'northeast', 'Fontsize', 18);
    end
    
    %% HORIZONTAL-FRAME VELOCITIES(subplotted)
    figure('Name','Horizontal Frame Velocities - All Flight','NumberTitle','off');
    
    % x axis
    subplot(3,1,1);
    Vhframe = [quatrotate(quatconj(Ya(:, 10:13)), Ya(:, 4:6)); Yf(Na + 1:end, 4:6)];
    plot(Tf, Vhframe(:, 1)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-x [m/s]');
    
    h = zeros(Np, 1);
    for i = 1:Np
        h(i) = plot(bound_value(i).t, bound_value(i).V(1), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    if not(settings.ballistic)
        legend(h(:), {'Apogee',  strcat('parachute ',  " " , string(2:Np), " ", 'opening')}, 'Location', 'southeast');
    else
        legend(h(:), {'Apogee'}, 'Location', 'southeast');
    end
    
    % y axis
    subplot(3,1,2)
    plot(Tf, Vhframe(:, 2)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-y [m/s]');
    
    for i = 1:Np
        h(i) = plot(bound_value(i).t, bound_value(i).V(2), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    % z axis
    subplot(3,1,3)
    plot(Tf, -Vhframe(:, 3)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-z [m/s]');
    
    for i = 1:Np
        h(i) = plot(bound_value(i).t, bound_value(i).V(2), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    %% ALTITUDE,MACH,VELOCITY,ACCELERATION(subplotted)
    figure('Name','Altitude, Mach, Velocity-Abs, Acceleration-Abs - ascent Phase','NumberTitle','off');
    subplot(2,3,1:3)
    plot(Ta, za), grid on, xlabel('time [s]'), ylabel('altitude [m]');
    
    subplot(2,3,4)
    plot(data_ascent.integration.t(1:end-1), data_ascent.interp.M(1:end-1)), grid on;
    xlabel('Time [s]'); ylabel('Mach M [/]')
    
    
    subplot(2,3,5)
    plot(Ta, abs_V), grid on;
    xlabel('time [s]'), ylabel('|V| [m/s]');
    
    
    subplot(2,3,6)
    plot(Ta, abs_A/9.80665), grid on;
    xlabel('time [s]'), ylabel('|A| [g]');
    
    
    %% TRAJECTORY PROJECTIONS(subplotted)
    figure('Name','Trajectory Projections - All Flight','NumberTitle','off');
    
    % xy frame
    subplot(1,3,1)
    plot(y, x), axis equal, hold on, grid on;
    xlabel('y, East [m]'), ylabel('x, North [m]');
    
    h = zeros(Np+2, 1);
    for i = 1:Np
        h(i) = plot(bound_value(i).X(1), bound_value(i).X(2), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    h(Np+1) = plot(0, 0, 'r.','markersize',14);
    h(Np+2) = plot(Yf(end,2), Yf(end,1), 'rx','markersize',7);
    
    if not(settings.ballistic)
            legend(h(:), {'Apogee', strcat('parachute ',  " " , string(2:Np), " ", 'opening'), 'Launch point',...
        'Landing point'}, 'Location', 'southeast');
    else
        legend(h(:), {'Apogee','Launch point', 'Landing point'}, 'Location', 'southeast');
    end
    
    % xz frame
    subplot(1,3,2)
    plot(x, z), hold on, grid on;
    xlabel('x, North [m]'), ylabel('z, Altitude [m]');
    
    for i = 1:Np
        plot(bound_value(i).X(2), bound_value(i).X(3), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    plot(0, 0, 'r.','markersize',14);
    plot(Yf(end, 1), -Yf(end, 3), 'rx','markersize',7);
    
    % yz frame
    subplot(1,3,3)
    plot(y, z), hold on, grid on;
    xlabel('y, East [m]'), ylabel('z, Altitude [m]');
    
    for i = 1:Np
        plot(bound_value(i).X(1), bound_value(i).X(3), 'o', 'MarkerSize',...
            10, 'MarkerFaceColor', Colors(:,i), 'MarkerEdgeColor', 'none');
    end
    
    plot(0, 0, 'r.','markersize',14);
    plot(Yf(end, 2), -Yf(end, 3), 'rx','markersize',7);
    
    %% EULERIAN ANGLES (subplotted)
    % pre-allocation
    pitch_angle = zeros(Na, 1); yaw_angle = zeros(Na, 1); roll_angle = zeros(Na, 1);
    
    % retriving the angle from the state
    for k = 2:Na
        pitch_angle(k) = pitch_angle(k-1) + (Ya(k, 8) + Ya(k-1, 8))/2*180/pi*(Ta(k) - Ta(k-1));
        yaw_angle(k) = yaw_angle(k-1) + (Ya(k, 9) + Ya(k-1, 9))/2*180/pi*(Ta(k) - Ta(k-1));
        roll_angle(k) = roll_angle(k-1) + (Ya(k, 7) + Ya(k-1, 7))/2*180/pi*(Ta(k) - Ta(k-1));
    end
    
    figure('Name', 'Eulerian Angles - ascent Phase', 'NumberTitle', 'off');
    subplot(3,1,1)
    plot(Ta, pitch_angle + settings.OMEGAmin*180/pi)
    grid on, xlabel('time [s]'), ylabel('pitch angle [deg]');
    
    subplot(3,1,2)
    plot(Ta, yaw_angle)
    grid on, xlabel('time [s]'), ylabel('yaw angle [deg]');
    
    subplot(3,1,3)
    plot(Ta, roll_angle)
    grid on, xlabel('time [s]'), ylabel('roll angle [deg]')
    
    
    
else   %%%% STOCHASTIC PLOTS (only if N>1)
    
    %% LANDING POINTS 2D
    % Position Scaled map in background
    figure('Name', 'Landing Points', 'NumberTitle','off')
    if settings.landing_map
        imshow(settings.map_file, 'YData', -settings.map_yaxis, 'XData', settings.map_xaxis);
        set(gca, 'YDir', 'normal'); % set y axis with ascending increasing values
    end
    axis on; hold on
    plot(LP(:, 2), LP(:, 1), '.r','MarkerSize', 11);
    title('Atterraggio con secondo drogue');
    xlabel('m')
    ylabel('m')
    axis image
    
    %% LANDING POINTS 3D
    
    figure('Name', 'Landing Points Map', 'NumberTitle', 'off')
    hold on; grid on;
    
    X_t = -8000:30:8000;
    Y_t = -8000:30:8000;
    L_X = length(X_t);
    L_Y = length(Y_t);
    Z_t = zeros(L_X, L_Y);
    for i = 1:L_X
        for j = 1:L_Y
            Z_t(j,i) = settings.funZ(X_t(i), Y_t(j));
        end
    end
    
    surf(X_t, -Y_t, -Z_t, 'EdgeColor', 'none'); colorbar; view(-90, 80);
    if not(settings.ballistic)
        h1 = plot3(LPin(:,1), -LPin(:,2), -LPin(:, 3),'o', 'MarkerEdgeColor', 'r', 'MarkerfaceColor', 'r');
        h2 = plot3(LPout(:, 1), -LPout(:, 2), -LPout(:, 3), 'o', 'MarkerEdgeColor', 'k', 'MarkerfaceColor', 'k');
        legend([h1, h2], 'safe landing points', 'not-safe landing points')
    else
        h1 = plot3(LP(:,1), -LP(:,2), -LP(:,3), 'o', 'MarkerEdgeColor', 'k', 'MarkerfaceColor', 'k');
        legend(h1,'ballistic landing points')
    end
    xlabel('North [m]'), ylabel('East [m]'), zlabel('Altitude [m]');
    
    %% LAST PARACHUTE OPENING POINTS
    if not(settings.ballistic)
        figure('Name', 'Last Drogue opening Points', 'NumberTitle', 'off')
        plot(0, 0, 'ro', 'MarkerSize', 20, 'MarkerFacecolor', 'b');
        hold on
        plot(RPin(:, 2),RPin(:, 1), 'o', 'MarkerEdgeColor', 'r', 'MarkerfaceColor', 'r');
        plot(RPout(:, 2),RPout(:, 1), 'o', 'MarkerEdgeColor', 'k', 'MarkerfaceColor', 'k');
        legend({'Launch Site', 'safe points', 'not-safe points'})
    end
    
    %% APOGEE POINTS
    figure('Name','Apogee Points','NumberTitle','off');
    plot(xapom, yapom, 'bs', 'MarkerSize', 20, 'MarkerFacecolor', 'b'), hold on;
    plot(X(:, 1), X(:, 2),'k+');
    legend('Mean Apogee Point', 'Apogee Points');
    title('Apogee Points'), xlabel('North [m]'), ylabel('East [m]');
    view(90, 270)
    axis equal
    
end
