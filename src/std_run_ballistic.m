function [ Tf,Yf, Ta,Ya ] = std_run_ballistic( settings )
%STD RUN - This function runs a standard (non-stochastic) simulation
% OUTPUTS
% Tf: Time steps
% Yf: Final State

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 31.XII.2014
% License:  2-clause BSD

% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

%Starting Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

%Starting State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];

% Wind Generation
[uw,vw,ww] = windgen(settings.wind.AzMin,settings.wind.AzMax,...
    settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
    settings.wind.MagMax);

%% ASCEND %%

[Ta,Ya] = ode45(@ascend,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww);

%% DESCEND %%

[Td,Yd] = ode45(@ballistic_descent,settings.ode.timedesc,Ya(end,:),settings.ode.optionsdesc,...
    settings,uw,vw,ww);

% %% DROGUE 1 %% 
% para = 1; %Flag for Drogue 1
% 
% %Initial Condition are the last from ascend (need to rotate because
% %velocities are outputted in body axes)
% X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
% [Td1,Yd1] = ode45(@descent_parachute,settings.ode.timedrg1,X0d1,...
%     settings.ode.optionsdrg1,settings,uw,vw,ww,para);
% 
% %% DROGUE 2 %% 
% para = 2; %Flag for Drogue 2
% 
% %Initial Condition are the last from drogue 1 descent
% X0d2 = Yd1(end,:);
% [Td2,Yd2] = ode45(@descent_parachute,settings.ode.timedrg2,X0d2,...
%     settings.ode.optionsdrg2,settings,uw,vw,ww,para);
% 
% %% MAIN %%
% para = 3; %Flag for Main (Rogall)
% 
% %Initial Condition are the last from drogue 2 descent
% X0m = Yd2(end,:);
% [Trog,Yrog] = ode45(@descent_parachute,settings.ode.timemain,X0m,...
%     settings.ode.optionsmain,settings,uw,vw,ww,para);

%% FINAL STATE ASSEMBLING %%

%Total State
Yf = [Ya;Yd];
%Total Time
Tf = [Ta; Ta(end)+Td];

%% PLOTTING THINGS %%

if settings.plot == 1
    
    set(0,'DefaultAxesFontSize',settings.DefaultFontSize,...
    'DefaultLineLineWidth',settings.DefaultLineWidth);

    %Interpolation for less points --> visualization
    Tinterp = linspace(Tf(1),Tf(end),settings.tSteps);
    [Tunique,ia,~]= unique(Tf);
    u = interp1(Tunique,Yf(ia,4),Tinterp);
    v = interp1(Tunique,Yf(ia,5),Tinterp);
    w = interp1(Tunique,Yf(ia,6),Tinterp);
    
    % VELOCITIES %
    Va=quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6)); %apogee
    
    figure;
    suptitle('Velocities Profiles on Time')
    subplot(3,1,1);
    plot(Tinterp,u,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Velocity-x [m/s]');
    h1=plot(Ta(end),Va(1),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    h2=plot(Ta(end)+Td(end),Yd(end,4),'ks','MarkerSize',8,...
        'MarkerFaceColor','k');
    legend([h1,h2],'Apogee/1st Drogue Deployment',...
        '2nd Drogue Deployment',...
        'Location','southeast');
    grid on
    
    subplot(3,1,2)
    plot(Tinterp,v,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Velocity-y [m/s]');
    h1=plot(Ta(end),Va(2),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    h2=plot(Ta(end)+Td(end),Yd(end,5),'ks','MarkerSize',8,...
        'MarkerFaceColor','k');
    legend([h1,h2],'Apogee/1st Drogue Deployment',...
        '2nd Drogue Deployment',...
        'Location','southeast');
    grid on
    
        
    subplot(3,1,3)
    plot(Tinterp,w,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Velocity-z [m/s]');
    h1=plot(Ta(end),Va(3),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    h2=plot(Ta(end)+Td(end),Yd(end,6),'ks','MarkerSize',8,...
        'MarkerFaceColor','k');
       legend([h1,h2],'Apogee/1st Drogue Deployment',...
        '2nd Drogue Deployment',...
        'Location','southeast');
    grid on
    
    %COMPLETE TRAJECTORY% 
    figure;
    h0=plot3(Yf(1,2),Yf(1,1),-Yf(1,3)+settings.z0,'k+','MarkerSize',10);
    hold on
    plot3(Yf(:,2),Yf(:,1),-Yf(:,3)+settings.z0,'k-','LineWidth',2)
    title('Complete Trajectory');
    xlabel('East [m]')
    ylabel('Noth [m]');
    zlabel('Altitude [m]');
    h1=plot3(Ya(end,2),Ya(end,1),-Ya(end,3)+settings.z0,'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    h2=plot3(Yd(end,2),Yd(end,1),-Yd(end,3)+settings.z0,'ks','MarkerSize',8,...
        'MarkerFaceColor','k');
    
    legend([h0,h1,h2],'Launch Pad','Apogee/1st Drogue Deployment',...
        '2nd Drogue Deployment','Main Parachute Deployment');
    grid on
    
    % PLANAR DISPLACEMENT %

    figure
    plot(Yf(:,2),Yf(:,1),'k-','LineWidth',2)
    hold on
    title('Planar Displacement');
    xlabel('East [m]')
    ylabel('North [m]');
    h1=plot(Ya(end,2),Ya(end,1),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    h2=plot(Yd(end,2),Yd(end,1),'ks','MarkerSize',8,...
        'MarkerFaceColor','k');
    
    legend([h1,h2],'Apogee/1st Drogue Deployment',...
        '2nd Drogue Deployment',...
        'Location','southeast');
    grid on
    
     
    
     
   %Interpolation for less points --> visualization
    Tinterp = linspace(Ta(1),Ta(end),settings.tSteps);
    p = interp1(Ta,Ya(:,7),Tinterp);
    q = interp1(Ta,Ya(:,8),Tinterp);
    r = interp1(Ta,Ya(:,9),Tinterp);
   
    % ANGULAR RATES %
    
    figure;
    suptitle('Angular rates on Time')
    subplot(3,1,1);
    plot(Tinterp,p,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Roll rate p [rad/s]');
    h1=plot(Ta(end),Ya(end,7),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    legend(h1,'Apogee/1st Drogue Deployment','Location','southeast');
    grid on
    
    subplot(3,1,2)
    plot(Tinterp,q,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Pitch rate q [rad/s]');
    h1=plot(Ta(end),Ya(end,8),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    legend(h1,'Apogee/1st Drogue Deployment','Location','northeast');
    grid on
    
        
    subplot(3,1,3)
    plot(Tinterp,r,'k-')
    hold on
    xlabel('Time[s]');
    ylabel('Yaw rate r [rad/s]');
    h1=plot(Ta(end),Ya(end,9),'ko','MarkerSize',8,...
        'MarkerFaceColor','k');
    legend(h1,'Apogee/1st Drogue Deployment','Location','southeast');
    grid on
    
end

%Resizing
h = get(0,'children');
scrsz = get(0,'ScreenSize');
for i=1:length(h)
  set(h(i),'OuterPosition',[0 0 scrsz(4) scrsz(4)])
  %saveas(h(i), ['figure' num2str(i)], 'fig');
end



end