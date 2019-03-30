clear all
close all
clc

% a = csvread('wind.txt', 0, 11, [0,11,908,11])
f = fopen('wind.txt','r');
a = textscan(f','%f %f %f %f %f %f %f %f %f %f %f %f %f');

speed = a{1,11}/3.6;
dir = a{1,13};

dir_reor = sort(dir);

[c,ia,ic] = unique(dir_reor);
    
ia(17) = 908;

for i = 1:length(ia)-1
    q(i) = ia(i+1) - ia(i);
end

[m, ind_max ] = max(q);

q = q/908*100;

date = 1:1:908;
date_date = zeros(908,6);
for i = 1:length(date)
    date_init = 6136;
    date_mjd2000 = date_init+date(i);
    date_date(i,1:6) = mjd20002date(date_mjd2000);
end

time = datenum(date_date);

mean(speed)
%%

labels = {'N','NNE','NE','ENE','E','ESE','SE','SSE','S',...
    'SSW','SW','WSW','W','WNW','NW','NNW'};
pie(q,labels)


% plot(dir,speed,'o');
% 
% hist(dir,speed)
% 
% 
figure()
a = polyfit(time,speed,3);
b = polyval(a,time(1):0.1:time(end));
plot(time(1):0.1:time(end),b)
hold on
plot(time,speed);
datetick('x',1,'keepticks')
% figure()
% plot(date,dir);
