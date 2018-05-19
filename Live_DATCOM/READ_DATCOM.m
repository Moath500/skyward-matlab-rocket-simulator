function coefficients=READ_DATCOM(alpha,Beta,Mach,Alt)

% Input of the timestep: This elements would not bedefined here, but they
% are an input for this script as a function, coming from the values at
% that time step.

alpha=1.5267894;
Beta=0.0585679;
Mach=0.4;
Alt=0;

%% Create alpha values for DATCOM computation:

Beta=floor(Beta*100)/100;

alphad=floor(alpha*10);
alphau=ceil(alpha*10);

alphav=[alphad-2, alphad-1, alphad, alphau, alphau+1, alphau+2];
alphav=alphav/10;

%% Create the for005.dat and for006.dat

% -------------------> Calling "Auto_Matrices.m"

% Auto_Matrices_sim

AERO=datcomimport('for006.dat');
AERO=AERO{1,1};

%% Take coeficients:

alphas=alphav(3:4);

coefficients.cl=interp1(alphas,AERO.cl(3:4),alpha);
coefficients.cd=interp1(alphas,AERO.cd(3:4),alpha);
coefficients.ca=interp1(alphas,AERO.ca(3:4),alpha);
coefficients.xcp=interp1(alphas,AERO.xcp(3:4),alpha);
coefficients.cma=interp1(alphas,AERO.cma(3:4),alpha);
coefficients.cyb=interp1(alphas,AERO.cyb(3:4),alpha);
coefficients.cnb=interp1(alphas,AERO.cnb(3:4),alpha);
coefficients.cna=interp1(alphas,AERO.cna(3:4),alpha);
coefficients.cmq=interp1(alphas,AERO.cmq(3:4),alpha);
coefficients.cnr=interp1(alphas,AERO.cnr(3:4),alpha);
coefficients.clp=interp1(alphas,AERO.clp(3:4),alpha);
coefficients.cmad=interp1(alphas,AERO.cmad(3:4),alpha);
coefficients.cnp=interp1(alphas,AERO.cnp(3:4),alpha);

end