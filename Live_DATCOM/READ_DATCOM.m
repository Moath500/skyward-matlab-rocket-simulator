% Author: Adriano Filippo Inno && Fernando Soler
% Skyward Experimental Rocketry | AFD Dept
% email: adriano.filippo.inno@skywarder.eu
% Release date: 20/05/2018

% This function given a single fligth condition computes the related
% cofficients using datcom and then datcomimport

% function Run_datcom needed


function coeff = READ_DATCOM(Alpha,Beta,Mach,Alt,Xcg)

%% SETTING DATCOM PRECISION

if abs(Alpha) < 2
    N_Alpha = 7;
    STEP_Alpha = 0.1;
elseif abs(Alpha) > 2 && abs(Alpha) < 5
    N_Alpha = 9;
    STEP_Alpha = 0.5;
else
    N_Alpha = 9;
    STEP_Alpha = 2;
end

%% CREATE FLIGHT CONDITION VALUES FOR DATCOM

Beta_datc = round(Beta,2);               % rounded to the 2nd decimal place
Mach_datc = round(Mach,2);               % rounded to the 2nd decimal place
Alt_datc = round(Alt);                   % rounded to the nearest integer
Alpha_centered = round(Alpha,2);         % rounded to the 2nd decimal place

Alpha_datc = Alpha_centered*ones(1,N_Alpha)+((-(N_Alpha-1)/2):((N_Alpha-1)/2))...
        *STEP_Alpha;  % vector of alphas needed in datcom
    

%% CREATE FOR005.DAT AND RUN DATCOM

Run_datcom(Alpha_datc,Beta_datc,Alt_datc,Mach_datc,Xcg)

AERO = datcomimport('for006.dat');
AERO = AERO{1,1};

%% INTERPOLATION

if Alpha > Alpha_datc((N_Alpha+1)/2)
    i = (N_Alpha+1)/2;
else
    i = (N_Alpha+1)/2-1;
end

Alpha_interp = Alpha_datc(i:(i+1));

coeff.cl   = interp1(Alpha_interp,AERO.cl(i:(i+1)),Alpha);
coeff.cd   = interp1(Alpha_interp,AERO.cd(i:(i+1)),Alpha);
coeff.ca   = interp1(Alpha_interp,AERO.ca(i:(i+1)),Alpha);
coeff.xcp  = interp1(Alpha_interp,AERO.xcp(i:(i+1)),Alpha);
coeff.cma  = interp1(Alpha_interp,AERO.cma(i:(i+1)),Alpha);
coeff.cyb  = interp1(Alpha_interp,AERO.cyb(i:(i+1)),Alpha);
coeff.cnb  = interp1(Alpha_interp,AERO.cnb(i:(i+1)),Alpha);
coeff.cna  = interp1(Alpha_interp,AERO.cna(i:(i+1)),Alpha);
coeff.cmq  = interp1(Alpha_interp,AERO.cmq(i:(i+1)),Alpha);
coeff.cnr  = interp1(Alpha_interp,AERO.cnr(i:(i+1)),Alpha);
coeff.clp  = interp1(Alpha_interp,AERO.clp(i:(i+1)),Alpha);
coeff.cmad = interp1(Alpha_interp,AERO.cmad(i:(i+1)),Alpha);
coeff.cnp  = interp1(Alpha_interp,AERO.cnp(i:(i+1)),Alpha);

end