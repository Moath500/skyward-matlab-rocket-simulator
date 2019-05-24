%{ 
Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
 Release date: 17/10/2018

This script computes allows to obtain the correlation between the wind
(magnitude and direction) and the rogallo deployment height in order to
keep the landing safe.

The Safe zone is here defined as an ellipse.
Remember to set random wind in config with 0-8 m/s , to fix the wind 
direction and to set upwind "true"

PARFOOR NEED DO BE STARTED BY HAND AND COMMENT IN STOCH_RUN.M THE COMMAND
"PARFOOR" near line 56.
%}



%%
clc
clear
close all

%%
hrog = 200:50:500;
p = zeros(length(hrog),2);

run('config_R2A_hermes.m');

settings.stoch.N = 100;


for k = 1:length(hrog)
    settings.zdrg2 = hrog(k);
    [LP,X,ApoTime,data_ascent,data_para] = stoch_run(settings);
    norm_wind = zeros(100,1);
    norm_R = zeros(100,1);
    for i = 1:100
        norm_wind(i) = norm(data_para{i}.wind(1).body_wind(1:3));
        norm_R(i) = norm(data_para{i}.state(2).Y(end,1:2));
    end
    
    [norm_wind,i] = sort(norm_wind);
    norm_R = norm_R(i);
    p(k,:) = polyfit(norm_wind,norm_R,1);
    
end

%%

figure()
imshow(settings.map_file, 'YData',-settings.map_yaxis, 'XData',settings.map_xaxis);
set(gca,'YDir','normal'); % set y axis with ascending increasing values
hold on
axis on

a = 1100;
b = 2800;
x0 = 0;
y0 = -300;
t = -pi:0.01:pi;
x = x0 + a*cos(t);
y = y0 + b*sin(t);
alpha = 10;

x_hat = cosd(alpha)*x - sind(alpha)*y;
y_hat = sind(alpha)*x + cosd(alpha)*y;
hold on
plot(x_hat,y_hat,'Linewidth',2)

xl = 10000*(rand(50,1) -0.5);
yl = 10000*(rand(50,1) -0.5);

% Find max R of ellipse (from launch point)
R = sqrt(x_hat.^2+y_hat.^2);
theta = atan2(y,x)*180/pi;

% Empirical law between landing distance and wind
% R_landing = c*w + r0

figure; hold on; grid on
w = zeros(length(hrog),length(t));
h = zeros(length(hrog),1);
for i = 1:length(hrog) 
    w(i,:) = (R-p(i,2))./p(i,1);

    h(i) = plot(mod(-(theta+alpha-90),360),w(i,:),'.');
end

legend(h(1:zeros(length(hrog))),strcat('hrog =',string(hrog)),'Location','Best','FontSize',18)
xlabel('Wind direction'); ylabel('Max wind'); xlim([0 360]);