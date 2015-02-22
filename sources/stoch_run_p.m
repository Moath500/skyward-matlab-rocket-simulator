function [ LP ,Z ] = stoch_run_p( settings )
%STD RUN - This function runs a stochastic simulation (parallel)
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
X0a = [X0;V0;W0;Q0;settings.m0;settings.Ixxf;settings.Iyyf;settings.Izzf];


%PreAllocation
LP = zeros(settings.stoch.N,3);
Z = zeros(settings.stoch.N,1);
ApoTime = zeros(settings.stoch.N,1);

parfor_progress(settings.stoch.N);
parpool;
parfor i=1:settings.stoch.N
    % Wind Generation
    [uw,vw,ww] = windgen(settings.wind.AzMin,settings.wind.AzMax,...
        settings.wind.ElMin,settings.wind.ElMax,settings.wind.MagMin,...
        settings.windd.MagMax);

    %% ASCEND %%

    [Ta,Ya] = ode45(@ascend,settings.ode.timeasc,X0a,settings.ode.optionsasc,...
        settings,uw,vw,ww);

    %% DROGUE %% 
    para = 1; %Flag for Drogue

    %Initial Condition are the last from ascend (need to rotate because
    %velocities are outputted in body axes)
    X0d = [Ya(end,1:3) quatrotate(quatconj(Ya(end,10:13)),Ya(end,4:6))];
    [Td,Yd] = ode45(@descent_parachute,settings.ode.timedrg,X0d,...
        settings.ode.optionsdrg,settings,uw,vw,ww,para);

    %% MAIN %%
    para = 2; %Flag for Main


    %Initial Condition are the last from drogue descent
    X0m = Yd(end,:);
    [Tm,Ym] = ode45(@descent_parachute,settings.ode.timemain,X0m,...
        settings.ode.optionsmain,settings,uw,vw,ww,para);

    %% FINAL STATE ASSEMBLING %%

    %Total State
    Yf = [Ya(:,1:3) quatrotate(quatconj(Ya(:,10:13)),Ya(:,4:6));Yd;Ym];
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
    
    
    plot(xm,ym,'bs','MarkerSize',20,'MarkerFacecolor','b');
    hold on
    
    %Point of launch
    plot(0,0,'ro','MarkerSize',20,'MarkerFacecolor','r');
   
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

delete(gcp);

%Resizing
h = get(0,'children');
scrsz = get(0,'ScreenSize');
for i=1:length(h)
  set(h(i),'OuterPosition',[0 0 scrsz(4) scrsz(4)])
  %saveas(h(i), ['figure' num2str(i)], 'fig');
end



end
