
load('ascend_plot.mat')

%% PLOTS

figure();
plot(ascend.t, ascend.T, '.'), title('Thrust vs time'), grid on;
xlabel('Time [s]'); ylabel('Thrust [N]');

figure();
plot(ascend.t, ascend.alpha*180/pi), title('alpha vs time'), grid on;
xlabel('Time [s]'); ylabel('alpha [deg]')

figure();
plot(ascend.t, ascend.beta*180/pi), title('beta vs time'), grid on;
xlabel('Time [s]'); ylabel('beta [deg]')

figure();
plot(ascend.t(1:end-1), ascend.M(1:end-1)), title('mach vs time'), grid on;
xlabel('Time [s]'); ylabel('Mach M [/]')

figure();
plot(ascend.t, ascend.CA), title('Aerodyn Coeff vs time'), grid on;
xlabel('Time [s]'); ylabel('Drag Coeff CD [/]')

figure();
plot(ascend.t, ascend.Drag), title('Drag vs time'), grid on;
xlabel('Time [s]'); ylabel('Drag D [N]')

figure();
plot(ascend.t, ascend.Forces), title('Axial Force vs time'), grid on;
xlabel('Time [s]'); ylabel('Axial force [N]')

figure();
plot(ascend.t, ascend.XCP,'.'), title('Stability margin vs time'), grid on;
xlabel('Time [s]'); ylabel('S.M.[/]')

