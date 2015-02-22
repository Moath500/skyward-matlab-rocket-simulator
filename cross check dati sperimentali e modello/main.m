% Launch 1
clear all
close all
clc

set(0,'DefaultLineLineWidth',3,'DefaultAxesFontSize',24',...
    'DefaultLineMarkerSize',12)

load altitude_1
load altitude_2
load simulation.mat

[x1,i1] = sort(altitude_1(:,1));
x1 = x1-x1(1);
[x2,i2] = sort(altitude_2(:,1));
x2 = x2-x2(1);
y1 = altitude_1(:,2);
y1 = y1(i1);
y2 = altitude_2(:,2);
y2 = y2(i2);


plot(x2,y2,'ks');
hold on
plot(T,-Y(:,3),'k-');
title('Experimental data - Model Cross Check')
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Flight Data','Simulation')

load velocity_2
[x1,i1] = sort(velocity_2(:,1));
x1 = x1-x1(1);
y1 = velocity_2(:,2);
y1 = y1(i1);

figure;
plot(x1,-y1,'ks');
hold on
plot(T,Y(:,6),'k-');
title('Experimental data - Model Cross Check')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('Flight Data','Simulation')


h = get(0,'children');
scrsz = get(0,'ScreenSize');
for i=1:length(h)
  set(h(i),'OuterPosition',[0 0 scrsz(4) scrsz(4)])
  %saveas(h(i), ['figure' num2str(i)], 'fig');
end
