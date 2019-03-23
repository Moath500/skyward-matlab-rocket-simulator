close all, clc

run('config_R2A_hermes');
% settings.map_file = 'map_roccaraso_2000.jpg'; % name of map for landing points
% settings.map_xaxis = [-2000 2000];  % limits for the data of the landing map
% settings.map_yaxis = [-2000 2000];

figure()
imshow(settings.map_file, 'YData',-settings.map_yaxis, 'XData',settings.map_xaxis);
set(gca,'YDir','normal'); % set y axis with ascending increasing values
hold on
axis on

a = 1000;
b = 3000;
x0 = 0;
y0= -500;
t = -pi:0.01:pi;
x = x0+a*cos(t);
y = y0+b*sin(t);
alpha = 12;

x_hat = cosd(alpha)*x-sind(alpha)*y;
y_hat = sind(alpha)*x+cosd(alpha)*y;
hold on
plot(x_hat,y_hat,'Linewidth',2)

xl = 10000*(rand(50,1) -0.5);
yl = 10000*(rand(50,1) -0.5);
% plot(xl,yl,'or','Linewidth',2);

% Find max R of ellipse (from launch point)
R = sqrt(x_hat.^2+y_hat.^2);
theta = atan2(y,x)*180/pi;

% Empirical law between landing distance and wind
% R_landing = c*w + r0
c = 363;
r0 = 2.16;
% Find max wind
w = (R - r0)./ c;

figure,
plot(mod(-(theta+alpha-90),360),w,'.'),
xlabel('Wind direction')
ylabel('Max wind');
xlim([0 360]);



