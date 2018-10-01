% Auto_Matrices

% Skyward Experimental Rocketry 

% Written by: | Mauro De Francesco | MSA Team Leader |
% mauro.defrancesco@skywarder.eu

% Approved by: | Mauro De Francesco | MSA Team Leader | 
% mauro.defrancesco@skywarder.eu 


% This program is used to determine the stability of a rocket through the
% change of some of the parameters of the fins.
% This require in the path for matlab, datcom.exe and datcom_parser.
% The program compiles a datcom input file, generate a datcom output file
% and converts it to a matrix.

%--------------------------------------------------------------------------

% Definition of the INPUT for the code:

% XcgF : Longitudinal position of Gravity Center in wet condition
% XcgE : Longitudinal position of Gravity Center in dry condition
% Sref : Reference area
% Lref : Longitudinal reference length
% Latref : Lateral reference length
% Tnose : Nose shape
% Lnose : Nose length
% Dnose : Nose diameter at base
% Lcenter : Centerbody length
% Dcenter : Centerbody diameter at base
% Dexit : Nozzle diameter for base drag calculation
% Base : Flag for base plume interaction
% Xle : Distance from missile nose to chord leading edge at each semi-span location
% Npanel : Numbers of panels
% Phif : Angle from each panel
% Ler : Leading edge radius at each span station
% Sta : Chord station used in measuring sweep
% Sspan : Semi-span location
% Chord : Panel chord at each semi-span location
% Zupper : Thickness to chord ratio of upper surface
% Lmaxu : Fraction of chord from leading edge to max thickness of upper surface
% Lflatu : Fraction of chord of constant thickness  section of upper surface

%--------------------------------------------------------------------------

tic
clear 
close all
clc

%% Flight Conditions

Mach=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
Alpha=[-20 -15 -10 -7.5 -5 -2.5 -1.5 -1 -0.5 0.0 0.5 1 1.5 2.5 5 7.5 10 15 20];
m=length(Mach);
n=length(Alpha);
Beta=0;
Alt=[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 8000 9000 10000];

%% Design Parameters

% Starting the loop for various dimension of the fins (change these if you
% need some other parameters)
radius= 0.087;
height= 0.15; %put the values for different height
LSspan=radius+height;
LChord=0.32; %put the values for different chords
d=0.05; %distance from end of rocket

%% Fixed Parameters

% Reference Quantities
Xcgf=2.624;
Xcge=2.429;
Sref=0.02378;
Lref=0.174;
Latref=0.174;

% Axisymmetric Body Geometry
Lnose=0.958;
Dnose=0.174;
Lcenter=3.442;
Dcenter=0.174;

% Define Fin Set n
Npanel=4;
Phif=[0 90 180 270];
Ler=2*0.0015;
Sta=0;
Sweep=64;

%% Creation of for005.dat

for L1=1:length(LSspan)
    for L2=1:length(LChord)
% Create folder for your case
        first_t=height(L1);
        second_t=LChord(L2);        
        Xle=4.4 - d - LChord(L2);
        Sspan=[0.087 LSspan(L1)];
        h_fin=Sspan(2)-Sspan(1);
        chord_out=h_fin*tand(20)+LChord(L2)-h_fin*tand(64);
        
        if chord_out>0
            newdir = sprintf('%g',first_t, -second_t, -d);
            mkdir(fullfile(newdir)); 
        end
        
        for Xcg=[Xcgf Xcge]

% Define Fin Set n (Variables)


if chord_out>0
    
Chord=[LChord(L2) chord_out];

Zupper=[2e-3/Chord(1) 2e-3/Chord(2)];
Lmaxu=[3e-2/Chord(1) 3e-2/Chord(2)];
Lflatu=[(Chord(1)-0.06)/Chord(1) (Chord(2)-0.06)/Chord(2)];


% Creating for005.dat file from previous data
fid = fopen(strcat(pwd,'\for005.dat'),'w+'); %w+ open or create file for reading and writing; discard existing contents
% Flight Conditions
fprintf(fid,'\n $FLTCON\r\n');
fprintf(fid, '  BETA=-5.,\r\n');
fprintf(fid, '  ALT=20*0.,\r\n');
fprintf(fid, '  NMACH=20.');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  MACH=');
fprintf(fid, '%.2f',Mach(1));
fprintf(fid, ',');
for i=2:10
    fprintf(fid, '%.1f',Mach(i));
    fprintf(fid, ',');
end
fprintf(fid, '%.2f',Mach(11));
fprintf(fid, ',');
fprintf(fid, '%.1f',Mach(12));
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid,  '  MACH(13)=');
fprintf(fid, '%.2f',Mach(13));
fprintf(fid, ',');
for i=14:20
    fprintf(fid, '%.1f',Mach(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  NALPHA=19.');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  ALPHA=');
for i=1:10
    fprintf(fid, '%.1f',Alpha(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  ALPHA(10)=');
for i=10:19
    fprintf(fid, '%.1f',Alpha(i));
    fprintf(fid, ',');
end
fprintf(fid, '$');

% Reference Quantities
fprintf(fid, '\r\n $REFQ\r\n');
fprintf(fid, '  XCG=');
fprintf(fid, '%.5f,\r\n', Xcg);
fprintf(fid, '  SREF=');
fprintf(fid, '%.5f,\r\n', Sref);
fprintf(fid, '  LREF=');
fprintf(fid, '%.3f,\r\n', Lref);
fprintf(fid, '  LATREF=');
fprintf(fid, '%.3f', Latref);
fprintf(fid, ',');
fprintf(fid, '$');

% Axisymmetric Body Geometry
fprintf(fid, '\r\n $AXIBOD\r\n');
fprintf(fid, '  TNOSE=KARMAN');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  LNOSE=');
fprintf(fid, '%.3f,\r\n', Lnose);
fprintf(fid, '  DNOSE=');
fprintf(fid, '%.3f,\r\n', Dnose);
fprintf(fid, '  LCENTR=');
fprintf(fid, '%.3f,\r\n', Lcenter);
fprintf(fid, '  DCENTR=');
fprintf(fid, '%.3f,\r\n', Dcenter);
fprintf(fid, '  DEXIT=0.161');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  BASE=.FALSE.,$');

% Finset
fprintf(fid, '\r\n $FINSET1\r\n');
fprintf(fid, '  XLE=');
fprintf(fid, '%.3f,\r\n', Xle);
fprintf(fid, '  NPANEL=');
fprintf(fid, '%.1f,\r\n', Npanel);
fprintf(fid, '  PHIF=');
for i=1:length(Phif)
fprintf(fid, '%.1f', Phif(i));
fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  LER=2*');
fprintf(fid, '%.4f,\r\n', Ler);
fprintf(fid, '  SWEEP=');
fprintf(fid, '%.1f,\r\n', Sweep);
fprintf(fid, '  STA=');
fprintf(fid, '%.1f,\r\n', Sta);
fprintf(fid, '  SSPAN=');
for i=1:length(Sspan)
fprintf(fid, '%.3f', Sspan(i));
fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  CHORD=');
for i=1:length(Chord)
fprintf(fid, '%.3f', Chord(i));
fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  SECTYP=HEX,');
fprintf(fid, '\r\n');
fprintf(fid, '  ZUPPER=');
for i=1:length(Zupper)
fprintf(fid, '%.4f', Zupper(i));
fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  LMAXU=');
for i=1:length(Lmaxu)
fprintf(fid, '%.4f', Lmaxu(i));
fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  LFLATU=');
for i=1:length(Lflatu)
fprintf(fid, '%.4f', Lflatu(i));
fprintf(fid, ',');
end
fprintf(fid, '$\r\n');

% Options
fprintf(fid, 'DERIV RAD\r\n');
fprintf(fid, 'DIM M\r\n');
fprintf(fid, 'DAMP\r\n');
fprintf(fid, 'SAVE\r\n');
fprintf(fid, 'NEXT CASE\r\n');

% Cases
for j=1:length(Alt)
    for k=1:length(Beta)
        if Beta(k)==-5 && Alt(j)==0
        else
        fprintf(fid,' $FLTCON\r\n');
        fprintf(fid,'  BETA=');
        fprintf(fid, '%.1f,\r\n', Beta(k));
        fprintf(fid,'  ALT=');
        fprintf(fid, '%.1f', Alt(j));
        fprintf(fid, '$\r\n');
        fprintf(fid, 'DERIV RAD\r\n');
        fprintf(fid, 'DIM M\r\n');
        fprintf(fid, 'DAMP\r\n');
        fprintf(fid, 'SAVE\r\n');
        fprintf(fid, 'NEXT CASE\r\n');
        end
    end
end
fclose(fid);
end        
end
 %% Creating .dat files+parsing

dos('datcoming')
pause(20)
dos('parsing')

value=0;
while value==0
    value=exist('for006.mat','file');
    pause(1);
end

%% Managing files

s=sprintf('%s', pwd,'\', newdir);

if Xcg==Xcgf
     movefile(strcat(pwd,'\for006.mat'),strcat(s,'\R2A_full.mat'));
     
else
    movefile(strcat(pwd,'\for006.mat'),strcat(s,'\R2A_empty.mat'));
end
    
end

        end

beep %gives an alarm when the program stopped working
toc