% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept
% email: adriano.filippo.inno@skywarder.eu
% Release date: 20/05/2018


% This function compiles the choesen datcom input(for005) file, for a single
% flight condition case and generates the datcom output(for006)


% This function require in the matlab path:

% - datcom.exe && datcoming_mac.bat && datcoming_mac_shortcut.lnk (MAC
% USER ONLY)

% where datcoming_mac.bat is a bash file that allows to start datcom and
% datcoming_mac_shortcut.lnk is a shortcut that HAS TO BE CREATED from the
% used virtual machine, and moved into the folder in order to let matlab  
% start the bat file.
% NOTE also that the virtual machine must be SWITCHED ON!

% the datcoming_mac.bat MUST be modified regardingly to the own matlab path




% Definition of the parameters:

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

function Run_datcom(Alpha,Beta,Alt,Mach,Xcg)

%% CONFIG

% Reference Quantities
Sref = 0.02378;
Lref = 0.174;
Latref = 0.174;

% Axisymmetric Body Geometry
Lnose = 0.958;
Dnose = 0.174;
Lcenter = 3.4552;
Dcenter = 0.174;
Xle = 3.8432;

% Define Fin Set n
radius = 0.087;
height = 0.15;
Sspan = [0.087 radius+height];
Chord = [0.32 0.067];
Npanel = 4;
Phif = [0 90 180 270];
Ler = 2*0.0015;
Sta = 0;
Sweep = 64;

%% FOR 005.DAT

Zupper = [2e-3/Chord(1) 2e-3/Chord(2)];
Lmaxu = [3e-2/Chord(1) 3e-2/Chord(2)];
Lflatu = [(Chord(1)-0.06)/Chord(1) (Chord(2)-0.06)/Chord(2)];
fid = fopen(strcat(pwd,'/for005.dat'),'w+'); %w+ open or create file for reading and writing; discard existing contents

% Flight Conditions
fprintf(fid,'\n $FLTCON\r\n');
fprintf(fid, '  BETA=');
fprintf(fid, '%.2f',Beta);
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  ALT=');
fprintf(fid, '%.2f',Alt);
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  NMACH=1.');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  MACH=');
fprintf(fid, '%.2f',Mach);
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  NALPHA=6.');
fprintf(fid, ',');
fprintf(fid, '\r\n');
fprintf(fid, '  ALPHA=');
for i = 1:6
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
for i = 1:length(Phif)
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
for i = 1:length(Sspan)
    fprintf(fid, '%.3f', Sspan(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  CHORD=');
for i = 1:length(Chord)
    fprintf(fid, '%.3f', Chord(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  SECTYP=HEX,');
fprintf(fid, '\r\n');
fprintf(fid, '  ZUPPER=');
for i = 1:length(Zupper)
    fprintf(fid, '%.4f', Zupper(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  LMAXU=');
for i = 1:length(Lmaxu)
    fprintf(fid, '%.4f', Lmaxu(i));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
fprintf(fid, '  LFLATU=');
for i = 1:length(Lflatu)
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
fclose(fid);


%% DATCOMING

if ismac   % mac procedure
    path = strcat(pwd,'/datcoming_mac_shortcut.lnk');
    command = strcat('open',{' '}, path);
    command = command{1};
    system(command);
else       % win procedure
    path = strcat(pwd,'/datcom.exe');
    command = strcat('strt',{' '}, path);
    command = command{1};
    system(command);
    
    value = 0;
    while value == 0
        value = exist('for006.dat','file');
        pause(0.05)
    end
    pause(0.5)
    
end







