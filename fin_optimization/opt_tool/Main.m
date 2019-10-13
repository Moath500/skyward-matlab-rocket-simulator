%%% Main %%%

clear 
close all
clc 

%% DATA
% principal data 
run config.m

% data of the analyis 
caso='matrici';

switch caso
    case 'z0'
        z0=[0 100 200 500];
        n=length(z0);
        var=z0;
    case 'ms'
        ms=4.3:0.1:7;
        m0=ms+settings.mp;
        n=length(ms);
        var=ms;
    case 'lrampa'
        lrampa=[3 4 5 6];
        n=length(lrampa);
        var=lrampa;
    case 'OMEGA'
        OMEGA=[85 86 87 89]*pi/180;
        n=length(OMEGA);
        var=OMEGA;
    case 'matrici'
        r_name={'R2A' 'R1X'};
        n=length(r_name);
        var=1:n;
end




%% RUN 

% preallocation 
grafico=true;

apogee=zeros(1,n);
max_a=zeros(1,n);
Vexit=zeros(1,n);


for i=1:n
    switch caso
        case 'z0'
            settings.z0=z0(i);
        case 'ms'
            settings.ms=ms(i);
            settings.m0=m0(i);
        case 'lrampa'
            settings.lrampa=lrampa(i);
        case 'OMEGA'
            settings.OMEGA=OMEGA(i);
        case 'matrici'
            filename=r_name{i};
            
            % Coefficients in full configuration
            filename_full = strcat(filename,'_full.mat');
            CoeffsF = load(filename_full,'Coeffs');
            settings.CoeffsF = CoeffsF.Coeffs;
            clear('CoeffsF');
            
            % Coefficients in empty configuration
            filename_empty = strcat(filename,'_empty.mat');
            CoeffsE = load(filename_empty,'Coeffs');
            settings.CoeffsE = CoeffsE.Coeffs;
            clear('CoeffsE');
            
            
            s = load(filename_full,'State');
            settings.Alphas = s.State.Alphas';
            settings.Betas = s.State.Betas';
            settings.Altitudes = s.State.Altitudes';
            settings.Machs = s.State.Machs';
            clear('s');
    end

    
    [apogee(i), max_a(i),Vexit(i),t,vect_XCP]=start_simulation(settings);
   
    % ottimizzazione 
    
    
    
    
    
    
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
plot(var,apogee,'o-')
xlabel(caso)
ylabel('apogee')

% plot max acceleration 
figure()
plot(var,max_a,'o-')
xlabel(caso)
ylabel('max |a|')







