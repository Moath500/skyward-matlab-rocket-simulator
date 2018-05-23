% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
% Release date: 13/01/2018

if settings.stoch.N == 1
    
    load('ascend_plot.mat') 
    
    %% SORTING VECTOR FOR ASCEND
    
    [ascend.t,I] = sort(ascend.t);
    ascend.T = ascend.T(I);
    ascend.Drag = ascend.Drag(I);
    ascend.Forces = ascend.Forces(I);
    ascend.alpha = ascend.alpha(I);
    ascend.beta = ascend.beta(I);
    ascend.CA = ascend.CA(I);
    ascend.XCP = ascend.XCP(I);
    
    
    %% ASCEND PLOTS
    
    if settings.only_XCP
        
        figure('Name','Stability Margin - Ascend Phase','NumberTitle','off');
        plot(ascend.t, -ascend.XCP,'.'), title('Stability margin vs time'), grid on;
        xlabel('Time [s]'); ylabel('S.M.[/]')
    else
        
        figure('Name','Stability Margin - Ascend Phase','NumberTitle','off');
        plot(ascend.t, -ascend.XCP,'.'), title('Stability margin vs time'), grid on;
        xlabel('Time [s]'); ylabel('S.M.[/]')
    
    figure('Name','Forces - Ascend Phase','NumberTitle','off');
    suptitle('Forces')
    subplot(3,1,1)
    plot(ascend.t, ascend.T, '.'), grid on;
    xlabel('Time [s]'); ylabel('Thrust [N]');
    
    subplot(3,1,2)
    plot(ascend.t, ascend.Drag),grid on;
    xlabel('Time [s]'); ylabel('Drag D [N]')
    
    subplot(3,1,3)
    plot(ascend.t, ascend.Forces), grid on;
    xlabel('Time [s]'); ylabel('Axial force [N]')
    
    
    figure('Name','Aerodynamics Angles - Ascend Ahase','NumberTitle','off');
    suptitle('Aerodynamics Angles')
    subplot(2,1,1)
    plot(ascend.t, ascend.alpha*180/pi), grid on;
    xlabel('Time [s]'); ylabel('alpha [deg]')
    
    subplot(2,1,2)
    plot(ascend.t, ascend.beta*180/pi), grid on;
    xlabel('Time [s]'); ylabel('beta [deg]')
    
    figure('Name','Drag Coefficient - Ascend Phase','NumberTitle','off');
    plot(ascend.t, ascend.CA), title('Drag Coefficient vs time'), grid on;
    xlabel('Time [s]'); ylabel('Drag Coeff CD [/]')
    
    
    %% 3D TRAJECTORY
    
    figure('Name','3D Trajectory - All Flight','NumberTitle','off');
    plot3(y,x,z), axis equal, hold on, grid on;
    title('Trajectory')
    xlabel('y, East [m]'), ylabel('x, North [m]'), zlabel('Altitude [m]')
    
    % concentric circles distance
    theta_plot = linspace(0,2*pi);
    R_plot = [1, 2, 3, 4, 5]*1000;
    for j = 1:length(R_plot)
        x_plot = R_plot(j)*cos(theta_plot');
        y_plot = R_plot(j)*sin(theta_plot');
        z_plot = zeros(length(theta_plot), 1);
        plot3(y_plot, x_plot, z_plot, '--r')
    end
    
    h1 = plot3(bound_value.Xd1(1),bound_value.Xd1(2),bound_value.Xd1(3),'ro',...
        'MarkerSize',7,'MarkerFaceColor','r');
    h4 = plot3(Y(end,2),Y(end,1),Y(end,3),'rx','markersize',7);
    h5 = plot3(0, 0, 0, '*');
    
    
    if settings.ballistic == false
        h2 = plot3(bound_value.Xd2(1),bound_value.Xd2(2),bound_value.Xd2(3),'rs',...
            'MarkerSize',7,'MarkerFaceColor','r');
        h3 = plot3(bound_value.Xrog(1),bound_value.Xrog(2),bound_value.Xrog(3),'rs',...
            'MarkerSize',7,'MarkerFaceColor','r');
        legend([h5,h1,h2,h3,h4],{'Launch point','Apogee/1st Drogue Deployment',...
            '2nd Drogue Deployment','Rogallo Deployment','Landing point'},'Location','northeast');
    else
        legend([h5,h1,h4],'Launch point','Apogee',...
            'Landing point','Location','northeast');
        if settings.sdf
            h2 = plot3(bound_value.Xd2(1),bound_value.Xd2(2),bound_value.Xd2(3),'rs',...
                'MarkerSize',7,'MarkerFaceColor','r');
            legend([h5,h1,h2,h4],'Launch point','Apogee','2nd Drogue Failure',...
                'Landing point','Location','northeast');
        end
    end
    
    
    
    %% STAGNATION TEMPERATURE
    
    figure('Name','Temperature Profile - All Flight','NumberTitle','off');
    plot(T, Tamb-273.15);
    hold on, grid on;
    plot(T, Ttot-273.15)
    title('Temperature Profile')
    xlabel('time [s]'), ylabel('Temperature [C]');
    legend('Surrounding', 'Total', 'location', 'best')
    
    %% HORIZONTAL-FRAME VELOCITIES(subplotted)
    
    figure('Name','Horizontal Frame Velocities - All Flight','NumberTitle','off');
    suptitle('Horizontal Frame Velocities')
    
    subplot(3,1,1);
    plot(T,Y(:,4)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-x [m/s]');
    
    h1 = plot(bound_value.td1,bound_value.Vd1(1),'ro','MarkerSize',7,...
        'MarkerFaceColor','r');
    
    if settings.ballistic == false
        h2 = plot(bound_value.td2,bound_value.Vd2(1),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h3 = plot(bound_value.trog,bound_value.Vrog(1),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3],'Apogee/1st Drogue Deployment','2nd Drogue Deployment',...
            'Rogallo Deployment','Location','southeast');
    else
        legend(h1,'Apogee','Location','southeast');
        if settings.sdf
            h2 = plot(bound_value.td2,bound_value.Vd2(1),'rs','MarkerSize',7,...
                'MarkerFaceColor','r');
            legend([h1,h2],'Apogee/1st Drogue Deployment','2nd Drogue Failure',...
                'Location','southeast');
        end
    end
    
    subplot(3,1,2)
    plot(T,Y(:,5)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-y [m/s]');
    
    h1 = plot(bound_value.td1,bound_value.Vd1(2),'ro','MarkerSize',7,...
        'MarkerFaceColor','r');
    
    if settings.ballistic == false
        h2 = plot(bound_value.td2,bound_value.Vd2(2),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h3 = plot(bound_value.trog,bound_value.Vrog(2),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3],'Apogee/1st Drogue Deployment','2nd Drogue Deployment',...
            'Rogallo Deployment','Location','southeast');
    else
        legend(h1,'Apogee','Location','southeast');
        if settings.sdf
            h2 = plot(bound_value.td2,bound_value.Vd2(2),'rs','MarkerSize',7,...
                'MarkerFaceColor','r');
            legend([h1,h2],'Apogee/1st Drogue Deployment','2nd Drogue Failure',...
                'Location','southeast');
        end
    end
    
    subplot(3,1,3)
    plot(T,-Y(:,6)), hold on, grid on, xlabel('Time[s]'), ylabel('Velocity-z [m/s]');
    
    h1 = plot(bound_value.td1,bound_value.Vd1(3),'ro','MarkerSize',7,...
        'MarkerFaceColor','r');
    
    if settings.ballistic == false
        
        h2 = plot(bound_value.td2,bound_value.Vd2(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h3 = plot(bound_value.trog,bound_value.Vrog(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3],'Apogee/1st Drogue Deployment','2nd Drogue Deployment',...
            'Rogallo Deployment','Location','southeast');
    else
        legend(h1,'Apogee','Location','southeast');
        if settings.sdf
            h2 = plot(bound_value.td2,bound_value.Vd2(3),'rs','MarkerSize',7,...
                'MarkerFaceColor','r');
            legend([h1,h2],'Apogee/1st Drogue Deployment','2nd Drogue Failure',...
                'Location','southeast');
        end
    end
    
    %% ALTITUDE,MACH,VELOCITY,ACCELERATION(subplotted)
    
    figure('Name','Altitude, Mach, Velocity-Abs, Acceleration-Abs - Ascend Phase','NumberTitle','off');
    suptitle('Altitude, Mach, Velocity-Abs, Acceleration-Abs')
    subplot(2,3,1:3)
    plot(Ta, z_a), grid on, xlabel('time [s]'), ylabel('altitude [m]');
    
    subplot(2,3,4)
    plot(ascend.t(1:end-1), ascend.M(1:end-1)), grid on;
    xlabel('Time [s]'); ylabel('Mach M [/]')
    
    
    subplot(2,3,5)
    plot(Ta, abs_Va), grid on;
    xlabel('time [s]'), ylabel('|V| [m/s]');
    
    
    subplot(2,3,6)
    plot(Ta, abs_Aa/9.80665), grid on;
    xlabel('time [s]'), ylabel('|A| [g]');
    
    
    %% TRAJECTORY PROJECTIONS(subplotted)
    
    
    figure('Name','Trajectory Projections - All Flight','NumberTitle','off');
    suptitle('Trajectory Projections')
    
    subplot(1,3,1)
    plot(y, x), axis equal, hold on, grid on;
    h1 = plot(Ya(end,2), Ya(end,1), 'r*','markersize',7);
    h2 = plot(0, 0, 'r.','markersize',14);
    h3 = plot(Y(end,2), Y(end,1), 'rx','markersize',7);
    xlabel('y, East [m]'), ylabel('x, North [m]');
    
    if settings.ballistic == false
        
        h4 = plot(bound_value.Xd2(1),bound_value.Xd2(2),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h5 = plot(bound_value.Xrog(1),bound_value.Xrog(2),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3,h4,h5],'Apogee/1st Drogue Deployment','Launch Point',...
            'Landing Point','2nd Drogue Deployment','Rogallo Deployment','Location','Best');
    else
        legend([h1,h2,h3],'Launch point','Apogee','Landing point',...
            'Location','northeast');
        if settings.sdf
            h4 = plot(bound_value.Xd2(1),bound_value.Xd2(2),'sr','MarkerSize',7,...
                'MarkerFaceColor','r');
            legend([h1,h2,h3,h4],'Apogee/1st Drogue Deployment','Launch Point',...
                'Landing Point','2nd Drogue Failure','Location','Best');
        end
    end
    
    
    subplot(1,3,2)
    plot(x, z), hold on, grid on;
    plot(Ya(end,1), -Ya(end,3), 'r*','markersize',7);
    plot(Y(end,1), -Y(end,3), 'rx','markersize',7);
    plot(0, 0, 'r.','markersize',14);
    xlabel('x, North [m]'), ylabel('z, Altitude [m]');
    
    if settings.ballistic == false
        
        h4 = plot(bound_value.Xd2(2),bound_value.Xd2(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h5 = plot(bound_value.Xrog(2),bound_value.Xrog(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3,h4,h5],'Apogee/1st Drogue Deployment','Launch Point',...
            'Landing Point','2nd Drogue Deployment','Rogallo Deployment','Location','Best');
    else

        if settings.sdf
            h4 = plot(bound_value.Xd2(2),bound_value.Xd2(3),'sr','MarkerSize',7,...
                'MarkerFaceColor','r');
        end
    end
    
    
    subplot(1,3,3)
    plot(y, z), hold on, grid on;
    plot(Ya(end,2), -Ya(end,3), 'r*','markersize',7);
    plot(Y(end,2), -Y(end,3), 'rx','markersize',7);
    plot(0, 0, 'r.','markersize',14);
    xlabel('y, East [m]'), ylabel('z, Altitude [m]');
    
    if settings.ballistic == false
        
        h4 = plot(bound_value.Xd2(1),bound_value.Xd2(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        h5 = plot(bound_value.Xrog(1),bound_value.Xrog(3),'sr','MarkerSize',7,...
            'MarkerFaceColor','r');
        legend([h1,h2,h3,h4,h5],'Apogee/1st Drogue Deployment','Launch Point',...
            'Landing Point','2nd Drogue Deployment','Rogallo Deployment','Location','Best');
    else
        
        if settings.sdf
            h4 = plot(bound_value.Xd2(1),bound_value.Xd2(3),'sr','MarkerSize',7,...
                'MarkerFaceColor','r');
        end
    end
    
    
    %% ANGULAR RATES (subplotted)
    
    figure('Name','Angular Rates - Ascend Phase','NumberTitle','off');
    suptitle('Angular Rates');
    subplot(3,1,1)
    plot(Ta, Ya(:, 8)*180/pi), grid on;
    xlabel('time [s]'), ylabel('pitch rate q [grad/s]')
    
    subplot(3,1,2)
    plot(Ta, Ya(:, 9)*180/pi), grid on;
    xlabel('time [s]'), ylabel('yaw rate r [grad/s]')
    
    subplot(3,1,3)
    plot(Ta, Ya(:, 7)*180/pi), grid on;
    xlabel('time [s]'), ylabel('roll rate p [grad/s]')
    
    
    %% EULERIAN ANGLES (subplotted)
    
    pitch_angle = zeros(length(Ta),1);
    yaw_angle = zeros(length(Ta),1);
    roll_angle = zeros(length(Ta),1);
    
    for k = 2:length(Ta)
        pitch_angle(k) = pitch_angle(k-1) + (Ya(k, 8) + Ya(k-1, 8))/2*180/pi*(T(k)-T(k-1));
        yaw_angle(k) = yaw_angle(k-1) + (Ya(k, 9) + Ya(k-1, 9))/2*180/pi*(T(k)-T(k-1));
        roll_angle(k) = roll_angle(k-1) + (Ya(k, 7) + Ya(k-1, 7))/2*180/pi*(T(k)-T(k-1));
    end
    
    figure('Name','Eulerian Angles - Ascend Phase','NumberTitle','off');
    suptitle('Eulerian Angles')
    subplot(3,1,1)
    plot(Ta, pitch_angle+settings.OMEGA*180/pi)
    grid on, xlabel('time [s]'), ylabel('pitch angle [deg]');
    
    subplot(3,1,2)
    plot(Ta, yaw_angle)
    grid on, xlabel('time [s]'), ylabel('yaw angle [deg]');
    
    subplot(3,1,3)
    plot(Ta, roll_angle)
    grid on, xlabel('time [s]'), ylabel('roll angle [deg]')
    end  
    
    
    %% STOCHASTIC PLOTS (only if N>1)
    
else

    
    % LANDING POINTS
    
    if not(settings.ao)
        figure('Name','Landing Points','NumberTitle','off');
        plot(xm,ym,'bs','MarkerSize',20,'MarkerFacecolor','b'), hold on;
        plot(LP(:,1),LP(:,2),'k+');
        plot(0,0,'ro','MarkerSize',20,'MarkerFacecolor','r');
        legend('Mean Landing Point','Landing Points','Launch Site');
        xlabel('North [m]'), ylabel('East [m]'),title('Landing Points');
    end
    
    
    % APOGEE POINTS
    
    figure('Name','Apogee Points','NumberTitle','off');
    plot(xapom,yapom,'bs','MarkerSize',20,'MarkerFacecolor','b'), hold on;
    plot(X(:,1),X(:,2),'k+');
    legend('Mean Apogee Point','Apogee Points');
    title('Apogee Points'),xlabel('North [m]'), ylabel('East [m]');
    view(90,270)
    axis equal
    
    
    % HISTOGRAM
    
    [f,x] = hist(X(:,3),10);
    figure('Name','Apogee Histogram','NumberTitle','off');
    bar(x,f/settings.stoch.N);
    title('Apogee Altitudes Distribution')
    xlabel('Apogee [m]')
    ylabel('n_i/n_{tot}')
    
end

clear variable
