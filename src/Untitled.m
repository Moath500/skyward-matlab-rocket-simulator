dataFile = zeros(length(data_ascent.integration.t)+length(data_bal.integration.t),4);
[~, ~, Press,~]= atmosisa([data_ascent.interp.alt'; data_bal.interp.alt' ]);

%% ACCELERATION IN BODY FRAME
t = data_ascent.state.T;
N = length(t);
ua = data_ascent.state.Y(:,4); % Velocity in INERTIAL frame
% Compute acceleration with finite differences
aaN = (ua(3:N)-ua(1:N-2))./(t(3:N)-t(1:N-2));
% Add derivatevs at boundaries
aaN = [ua(2)/t(2); aaN; (ua(end)-ua(end-1))/(t(end)-t(end-1))];
% Add static acceleration (in BODY frame)
as = quatrotate(data_ascent.state.Y(:,10:13),[0 0 -9.806]);
% Recall: a_B = a_N - omega ^ V 
cor = cross(data_ascent.state.Y(:,7:9),data_ascent.state.Y(:,4:6));
aa = aaN - cor(:,1) + as(:,1);

ub = data_bal.state.Y(:,4);
N = length(ub);
t = data_bal.state.T;
abN = (ub(3:N)-ub(1:N-2))./(t(3:N)-t(1:N-2));
abN = [ub(2)/t(2); abN; (ub(end)-ub(end-1))/(t(end)-t(end-1))];
as = quatrotate(data_bal.state.Y(:,10:13),[0 0 -9.806]);
cor = cross(data_bal.state.Y(:,7:9),data_bal.state.Y(:,4:6));
ab = abN - cor(:,1) + as(:,1);

ax = [aa; ab]./9.80665;
dataFile(:,1) = T;
dataFile(:,2) = [-data_ascent.state.Y(:,3); -data_bal.state.Y(:,3)];
dataFile(:,3) = Press';
dataFile(:,4) = ax;


fileID = fopen('ballistic2.txt','w');
formatSpec = '%4.2f %4.2f %6.2f %4.2f \r\n';
fprintf(fileID,'%4s %4s %6s %6s \n','t[s]','h[m]','P[Pa]','ax[g]');
fprintf(fileID,formatSpec,dataFile');
fclose(fileID);


%% PLOT DESCENT
% Acceleration
% figure();
% subplot(1,2,1);
% plot(data_bal.state.T(1:end-1), abN(1:end-1)./9.80665), grid on, 
% xlabel('time [s]'), ylabel('ax Inertial [g]');
% subplot(1,2,2);
% plot(data_bal.state.T(1:end-1), ab(1:end-1)./9.80665), grid on, 
% xlabel('time [s]'), ylabel('ax Body [g]');

figure();
subplot(1,2,1);
plot(T, ax), grid on, 
xlabel('time [s]'), ylabel('ax [g]');
subplot(1,2,2);
plot(T, [aaN ;abN]./9.806), grid on, 
xlabel('time [s]'), ylabel('ax [g]');


% Inertial velocity
data_bal.state.Y(:,4:6) = quatrotate(quatconj(data_bal.state.Y(:,10:13)),data_bal.state.Y(:,4:6));
figure();
subplot(1,3,1)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,4)), grid on, 
xlabel('time [s]'), ylabel('u [m/s]');

subplot(1,3,2)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,5)), grid on, 
xlabel('time [s]'), ylabel('v [m/s]');

subplot(1,3,3)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,6)), grid on, 
xlabel('time [s]'), ylabel('w [m/s]');

% Angular velocity

figure();
subplot(1,3,1)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,7).*180/pi), grid on, 
xlabel('time [s]'), ylabel('p [deg/s]');

subplot(1,3,2)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,8).*180/pi), grid on, 
xlabel('time [s]'), ylabel('q [deg/s]');

subplot(1,3,3)
plot(data_bal.state.T(1:end-1), data_bal.state.Y(1:end-1,9).*180/pi), grid on, 
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
subplot(1,2,1)
plot(data_bal.state.T(1:end-1),data_bal.interp.alpha(1:end-1).*180/pi);
xlabel('time [s]'); ylabel('alpha [deg]');

subplot(1,2,2)
plot(data_bal.state.T(1:end-1),data_bal.interp.beta(1:end-1).*180/pi);
xlabel('time [s]'); ylabel('beta [deg]');

% Angle
pitch_angle = zeros(length(data_bal.state.T),1);
yaw_angle = zeros(length(data_bal.state.T),1);
roll_angle = zeros(length(data_bal.state.T),1);

for k = 2:length(data_bal.state.T)
    pitch_angle(k) = pitch_angle(k-1) + (data_bal.state.Y(k, 8) + ...
        data_bal.state.Y(k-1, 8))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
    yaw_angle(k) = yaw_angle(k-1) + (data_bal.state.Y(k, 9) + ...
        data_bal.state.Y(k-1, 9))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
    roll_angle(k) = yaw_angle(k-1) + (data_bal.state.Y(k, 7) + ...
        data_bal.state.Y(k-1, 7))/2*180/pi*(data_bal.state.T(k)-data_bal.state.T(k-1));
end

figure();
subplot(3,1,1)
plot(data_bal.state.T, pitch_angle +3.45)
grid on, xlabel('time [s]'), ylabel('pitch angle [deg]');

subplot(3,1,2)
plot(data_bal.state.T, yaw_angle + 4.34)
grid on, xlabel('time [s]'), ylabel('yaw angle [deg]');

subplot(3,1,3)
plot(data_bal.state.T, roll_angle + 0.41)
grid on, xlabel('time [s]'), ylabel('roll angle [deg]')
