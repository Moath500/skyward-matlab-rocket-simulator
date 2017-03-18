function [ LP ,Z ] = stoch_run( settings )
%STD RUN - This function runs a stochastic simulation
% OUTPUTS
% LP: Landing Points
% Z: Apogee Altitudes

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 29.V.2014
% License:  2-clause BSD

%Starting Attitude
Q0 = angle2quat(settings.PHI,settings.OMEGA,0*pi/180,'ZYX')';

%Starting State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
m_tot = settings.m0 + settings.para1.mass + settings.para2.mass + settings.para3.mass; %Overall Mass
X0a = [X0;V0;W0;Q0;m_tot;settings.Ixxf;settings.Iyyf;settings.Izzf];

%PreAllocation
LP = zeros(settings.stoch.N,3);
Z = zeros(settings.stoch.N,1);
ApoTime = zeros(settings.stoch.N,1);


parfor_progress(settings.stoch.N);
for i=1:settings.stoch.N
    %% Wind Generation
if settings.hwm
    uw = 0;
    vw = 0;
    ww = 0;
else
    [uw,vw,ww] = windgen(settings.wind.AzMin,settings.wind.AzMax,...
    settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
    settings.wind.MagMax);
end
%% ASCEND %%

[Ta,Ya] = ode45(@ascend,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
    settings,uw,vw,ww);

%% DROGUE 1 %% 
para = 1; %Flag for Drogue 1

%Initial Condition are the last from ascend (need to rotate because
%velocities are outputted in body axes)
X0d1 = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
[Td1,Yd1] = ode45(@descent_parachute,settings.ode.timedrg1,X0d1,...
    settings.ode.optionsdrg1,settings,uw,vw,ww,para);

%% DROGUE 2 %% 
para = 2; %Flag for Drogue 2

%Initial Condition are the last from drogue 1 descent
X0d2 = Yd1(end,:);
[Td2,Yd2] = ode45(@descent_parachute,settings.ode.timedrg2,X0d2,...
    settings.ode.optionsdrg2,settings,uw,vw,ww,para);

%% MAIN %%
para = 3; %Flag for Main (Rogall)

%Initial Condition are the last from drogue 2 descent
X0m = Yd2(end,:);
[Trog,Yrog] = ode45(@descent_parachute,settings.ode.timemain,X0m,...
    settings.ode.optionsmain,settings,uw,vw,ww,para);

    %% FINAL STATE ASSEMBLING %%

    %Total State
    Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd1; Yd2;Yrog];
    %Total Time
    %Tf = [Ta; Ta(end)+Td;Ta(end)+Td(end)+Tm];

    LP(i,:) = Yf(end,1:3);
    Z(i) = -Ya(end,3);
    ApoTime(i) = Ta(end);

    parfor_progress;

end

%Checking bad simulations
if numel(LP(LP(:,3)<-10,:))
    fprintf(['Some landing points might be incorrect' ...
        'please check parameters!\n']);
end

% Writing Things

%Mean Landing Point
xm = mean(LP(:,1));
ym = mean(LP(:,2));

%Mean Apogee Time
ApoTimem = mean(ApoTime);

%Std. Deviation Apogee Time
ApoTimestd = std(ApoTime);

%Mean Altitude
Zm = mean(Z);

%Std. Deviation Altitude
Zstd = std(Z);


% Printing to screen
text =['Mean Landing Point:X:%3.3f m, Y:%3.3f m\n',...
    'Mean Altitude: %3.3f m || STD: %3.3f m\n',...
    'Mean Apogee Time: %3.3f s || STD: %3.3f s\n'];
fprintf(text,xm,ym,Zm,Zstd,ApoTimem,ApoTimestd);

if settings.plot == 1
    %% PLOTTING THINGS
    
    
    plot(xm,ym,'ks','MarkerSize',20);
    hold on
    
    %Point of launch
    plot(0,0,'ro','MarkerSize',20);
   
    %All the landing points
    plot(LP(:,1),LP(:,2),'k+');
     
    title('Landing Points');
    xlabel('X [m]');
    ylabel('Y [m]');
    legend('Mean Landing Point','Launch Site','Landing Points');
    
    
    %Histogram
    [f,x] = hist(Z,10);
    figure;
    bar(x,f/settings.stoch.N);
    title('Apogee Altitudes Distribution')
    xlabel('Apogee [m]')
    ylabel('n_i/n_{tot}')

    


end

%Resizing
h = get(0,'children');
scrsz = get(0,'ScreenSize');
for i=1:length(h)
  set(h(i),'OuterPosition',[0 0 scrsz(4) scrsz(4)])
  %saveas(h(i), ['figure' num2str(i)], 'fig');
end




end
