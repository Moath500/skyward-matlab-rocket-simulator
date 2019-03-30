clear all
close all
clc

%% Data settings

% Data from .txt file retrieved from
% https://www.caputfrigoris.it/rete-meteo/roccaraso/daily.htm
fid = fopen('wind.txt','r');
data = textscan(fid','%f %f %f %f %f %f %f %f %f %f %f %f %f');

% Extracting wind data
speed = data{1,11}/3.6; % Speed [m/s]
dir = data{1,13}; % Direction [deg]

% Setting data for the date
date = 1:1:908;
date_date = zeros(908,6);
for i = 1:length(date)
    date_init = 6136;
    date_mjd2000 = date_init+date(i);
    date_date(i,1:6) = mjd20002date(date_mjd2000);
end
time = datenum(date_date);

%% Pie chart for the directions

dir_reor = sort(dir); % Sorts the directions
[c,ia,ic] = unique(dir_reor); 
ia(17) = 908;
q = zeros(length(ia)-1,1);
for i = 1:length(ia)-1
    q(i) = ia(i+1) - ia(i);
end

[m, ind_max ] = max(q);
q = q/908*100;

labels = {'N','NNE','NE','ENE','E','ESE','SE','SSE','S',...
    'SSW','SW','WSW','W','WNW','NW','NNW'};
explode = ones(length(q),1);
H = pie(q,explode);
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),labels(:))
title('Wind direction probability','FontSize',15,'FontWeight','bold')
legend('N: Nord [0 Deg]','NNE: Nord-Nord-Est [0-45 deg]','NE: Nord-Est [45 deg]',...
    'ENE: Est-Nord-Est [45-90 deg]','E: Est [90 deg]','ESE: Est-Sud-Est [90-135 deg]',...
    'SE: Sud-Est [135 deg]','SSE: Sud-Sud-Est [135-180 deg]',...
    'S: Sud [180 deg]','SSW: Sud-Sud-Ovest [180-225 deg]',...
    'SW: Sud-Ovest [225 deg]','WSW: Ovest-Sud-Ovest [225-270 deg]',...
    'W: Ovest [270 deg]','WNW: Ovest-Nord-Ovest [270-315 deg]',...
    'NW: Nord-Ovest [315 deg]','NNW: Nord-Nord-Ovest [315-360 deg]',...
    'location','best','FontSize',12)

%% Speed data

figure()
hold on
grid minor
mean_speed = mean(speed);

% Plot speed
plot(time,speed,'k');

% Plot mean of the speed
line([time(1),time(end)],[mean_speed,mean_speed],... 
    'Color','red','LineStyle','--','LineWidth',2);

% Plot polynomial fitting of speed
data = polyfit(time,speed,3);
b = polyval(data,time(1):0.1:time(end));
plot(time(1):0.1:time(end),b,'b','LineWidth',2)

datetick('x',1,'keepticks')
xlim([time(1) time(end)])
xlabel('Date','FontSize',15,'FontWeight','bold')
ylabel('Speed [m/s]','FontSize',15,'FontWeight','bold')
title('Hystory of wind speed in Roccaraso','FontSize',15,'FontWeight','bold')
legend('Day-to-day speed','Mean speed','Polynomial fit','FontSize',12)