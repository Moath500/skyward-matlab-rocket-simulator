%%% Main %%%

clear 
close all
clc 

%% DATA
% principal data 
run config.m

% data of the analyis 

ms=4.3:0.1:7;
m0=ms+settings.mp;
n=length(ms);



%% RUN 

% preallocation 
grafico=false;

apogee=zeros(1,n);
max_a=zeros(1,n);
Vexit=zeros(1,n);


for i=1:n
 
    settings.ms=ms(i);
    settings.m0=m0(i);

    [apogee(i), max_a(i),Vexit(i),t,vect_XCP]=start_simulation(settings);
   
   
    % grafici stabilità 
    
    if grafico 
        plot(t,vect_XCP,'.')
        xlabel('time [s]')
        ylabel('S / M')
        hold on 
    end
    
    
end


%% PLOT 

% plot apogee
figure()
plot(ms,apogee,'o-')
xlabel('structural mass [kg]')
ylabel('apogee [m]')

% plot max acceleration 
figure()
plot(ms,max_a,'o-')
xlabel('structural mass [kg]')
ylabel('max |a| [g]')







