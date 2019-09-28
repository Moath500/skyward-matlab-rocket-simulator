dataFile = zeros(length(data_ascent.integration.t)+length(data_bal.integration.t),8);
[~, ~, Press,~]= atmosisa([data_ascent.interp.alt'; data_bal.interp.alt' ]);

%% ACCELERATION IN BODY FRAME
t = data_ascent.state.T;
N = length(t);
ua = data_ascent.state.Y(:,4); % Velocity in body frame
% Compute acceleration with finite differences
aaN = (ua(3:N)-ua(1:N-2))./(t(3:N)-t(1:N-2));
% Add derivatevs at boundaries
aaN = [ua(2)/t(2); aaN; (ua(end)-ua(end-1))/(t(end)-t(end-1))];
% Add static acceleration (in BODY frame)
as = quatrotate(data_ascent.state.Y(:,10:13),[0 0 -9.806]);
aa = aaN + as(:,1);
aa_m = aa - cross(data_ascent.state.Y(:,7:9),data_ascent.state.Y(:,4:6));

ub = data_bal.state.Y(:,4);
N = length(ub);
t = data_bal.state.T;
abN = (ub(3:N)-ub(1:N-2))./(t(3:N)-t(1:N-2));
abN = [ub(2)/t(2); abN; (ub(end)-ub(end-1))/(t(end)-t(end-1))];
as = quatrotate(data_bal.state.Y(:,10:13),[0 0 -9.806]);
ab = abN + as(:,1);
ab_m = ab - cross(data_bal.state.Y(:,7:9),data_bal.state.Y(:,4:6));

ax = [aa; ab]./9.80665;

%% ANGLES DURING ASCENT
pitch_angle = zeros(length(data_ascent.state.T),1);
yaw_angle = zeros(length(data_ascent.state.T),1);
roll_angle = zeros(length(data_ascent.state.T),1);
pitch_angle(1) = 80;

for k = 2:length(data_ascent.state.T)
    pitch_angle(k) = pitch_angle(k-1) + (data_ascent.state.Y(k, 8) + ...
        data_ascent.state.Y(k-1, 8))/2*180/pi*(data_ascent.state.T(k)-data_ascent.state.T(k-1));
    yaw_angle(k) = yaw_angle(k-1) + (data_ascent.state.Y(k, 9) + ...
        data_ascent.state.Y(k-1, 9))/2*180/pi*(data_ascent.state.T(k)-data_ascent.state.T(k-1));
    roll_angle(k) = roll_angle(k-1) + (data_ascent.state.Y(k, 7) + ...
        data_ascent.state.Y(k-1, 7))/2*180/pi*(data_ascent.state.T(k)-data_ascent.state.T(k-1));
end
pitch0 = pitch_angle(end);
yaw0 = yaw_angle(end);
roll0 = roll_angle(end);

%% PLOT 
% Acceleration
% figure();
% subplot(1,2,1);
% plot(data_bal.state.T(1:end-1), abN(1:end-1)./9.80665), grid on, 
% xlabel('time [s]'), ylabel('ax Inertial [g]');
% subplot(1,2,2);
% plot(data_bal.state.T(1:end-1), ab(1:end-1)./9.80665), grid on, 
% xlabel('time [s]'), ylabel('ax Body [g]');

figure(); hold on, grid on
% plot(T, [aaN ;abN]./9.806),
plot(T, ax)
ax_m = [aa_m ; ab_m]./9.806;
% plot(T,ax_m(:,1)),
xlabel('time [s]'), ylabel('ax [g]');
% legend('FD','1g def','1g def + centr');


% % Velocity
% % Inertial
% data_bal.state.Y(:,4:6) = quatrotate(quatconj(data_bal.state.Y(:,10:13)),data_bal.state.Y(:,4:6));
% figure();
% subplot(1,3,1)
% plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,4)), grid on, 
% xlabel('time [s]'), ylabel('u [m/s]');
% 
% subplot(1,3,2)
% plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,5)), grid on, 
% xlabel('time [s]'), ylabel('v [m/s]');
% 
% subplot(1,3,3)
% plot(data_bal.state.T(1:end-1), -data_bal.state.Y(1:end-1,6)), grid on, 
% xlabel('time [s]'), ylabel('w [m/s]');

% Angular velocity

figure();
subplot(1,3,1), hold on
plot(T(1:end-1), [data_ascent.state.Y(:,7); data_bal.state.Y(1:end-1,7)].*180/pi), grid on, 
h1 = plot(data_ascent.state.T(end),data_ascent.state.Y(end,7)*180/pi,'or');
xlabel('time [s]'), ylabel('p [deg/s]');
legend(h1,'Apogee');

subplot(1,3,2), hold 
plot(T(1:end-1),[data_ascent.state.Y(:,8); data_bal.state.Y(1:end-1,8)].*180/pi), grid on,
h1 = plot(data_ascent.state.T(end),data_ascent.state.Y(end,8)*180/pi,'or');
legend(h1,'Apogee');
xlabel('time [s]'), ylabel('q [deg/s]');

subplot(1,3,3), hold on
plot(T(1:end-1), [data_ascent.state.Y(:,9); data_bal.state.Y(1:end-1,9)].*180/pi), grid on, 
h1 = plot(data_ascent.state.T(end),data_ascent.state.Y(end,9)*180/pi,'or');
legend(h1,'Apogee');
xlabel('time [s]'), ylabel('r [deg/s]');

% Body forces

% figure();
% subplot(1,3,1)
% plot(data_bal.state.T(1:end-1), data_bal.forces.AeroDyn_Forces(1,1:end-1)), grid on, 
% xlabel('time [s]'), ylabel('X-body force [N]');
% 
% subplot(1,3,2)
% plot(data_bal.state.T(1:end-1), data_bal.forces.AeroDyn_Forces(2,1:end-1)), grid on, 
% xlabel('time [s]'), ylabel('Y-body force [N]');
% 
% subplot(1,3,3)
% plot(data_bal.state.T(1:end-1), data_bal.forces.AeroDyn_Forces(3,1:end-1)), grid on, 
% xlabel('time [s]'), ylabel('Z body force [N]');

% Aerodynamic angles
figure();
subplot(1,2,1), hold on
plot(T(1:end-1),[data_ascent.interp.alpha data_bal.interp.alpha(1:end-1)].*180/pi);
h1 = plot(data_ascent.state.T(end),data_ascent.interp.alpha(end).*180/pi,'or');
legend(h1,'Apogee');
xlabel('time [s]'); ylabel('alpha [deg]');

subplot(1,2,2), hold on
plot(T(1:end-1),[data_ascent.interp.beta data_bal.interp.beta(1:end-1)].*180/pi);
h1 = plot(data_ascent.state.T(end),data_ascent.interp.beta(end).*180/pi,'or');
legend(h1,'Apogee');
xlabel('time [s]'); ylabel('beta [deg]');

% Angle
pitch_angle_d = zeros(length(data_bal.state.T),1);
yaw_angle_d = zeros(length(data_bal.state.T),1);
roll_angle_d = zeros(length(data_bal.state.T),1);

pitch_angle_d(1) = pitch0;
yaw_angle_d(1) = yaw0;
roll_angle_d(1) = roll0;

for k = 2:length(data_bal.state.T)
    pitch_angle_d(k) = pitch_angle_d(k-1) + (data_bal.state.Y(k, 8) + ...
        data_bal.state.Y(k-1, 8))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
    yaw_angle_d(k) = yaw_angle_d(k-1) + (data_bal.state.Y(k, 9) + ...
        data_bal.state.Y(k-1, 9))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
    roll_angle_d(k) = roll_angle_d(k-1) + (data_bal.state.Y(k, 7) + ...
        data_bal.state.Y(k-1, 7))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
end

figure();
subplot(3,1,1), hold on
plot(T, [pitch_angle; pitch_angle_d]);
h1 = plot(data_ascent.state.T(end),pitch_angle(end),'or');
grid on, xlabel('time [s]'), ylabel('pitch angle [deg]');
legend(h1, 'Apogee');

subplot(3,1,2), hold on
plot(T, [yaw_angle; yaw_angle_d]);
h1 = plot(data_ascent.state.T(end),yaw_angle(end),'or');
grid on, xlabel('time [s]'), ylabel('yaw angle [deg]');
legend(h1, 'Apogee');

subplot(3,1,3), hold on
plot(T, [roll_angle; roll_angle_d]);
h1 = plot(data_ascent.state.T(end),roll_angle(end),'or');
grid on, xlabel('time [s]'), ylabel('roll angle [deg]')
legend(h1, 'Apogee');


%% SAVE TO LOG FILE
dataFile(:,1) = T;
dataFile(:,2) = [-data_ascent.state.Y(:,3); -data_bal.state.Y(:,3)];
dataFile(:,3) = Press';
dataFile(:,4) = ax;
dataFile(:,5) = [data_ascent.interp.alpha data_bal.interp.alpha].*180/pi;
dataFile(:,6) = [data_ascent.interp.beta data_bal.interp.beta].*180/pi;
dataFile(:,7) = [pitch_angle; pitch_angle_d];
dataFile(:,8) = [yaw_angle; yaw_angle_d];
dataFile(:,9) = [roll_angle; roll_angle_d];


fileID = fopen('ballistic2.txt','w');
formatSpec = '%4.3f %4.3f %8.2f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\r\n';
fprintf(fileID,'%4s %4s %6s %6s %5s %5s %5s %5s %5s \r\n','t[s]','h[m]',...
    'P[Pa]','ax[g]','alfa[deg]','beta[deg]','pitch[deg]','yaw[deg]','rol[deg]');
fprintf(fileID,formatSpec,dataFile');
fclose(fileID);