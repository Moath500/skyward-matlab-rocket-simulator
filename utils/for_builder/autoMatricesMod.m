function [Coeffs, State] = autoMatricesMod(s)
%% States
% State values in which the aerodynamic coefficients will be computed
Mach = s.Mach;
Alpha = s.Alpha;
Alt = s.Alt;
Beta = s.Beta;
Nm = length(Mach);
Na = length(Alpha);
Nb = length(Beta);
Nalt = length(Alt);

shape = s.shape;

xcg = s.xcg;

C1 = s.C1;
C2 = s.C2;
H = s.H;


%% Fixed Parameters
D = 0.09;                                           % [m] rocket diameter
r = D/2;                                            % [m] rocket radius
S = r^2*pi;                                         % [m^2] cross section                         
Lnose = 0.3;                                        % [m] nose length
Lcenter = 1.72;                                     % [m] Lcenter : Centerbody length
Npanel = 4;                                         % [m] number of fins
Phif = [0 90 180 270];                              % [deg] Angle of each panel
Ler = 0.003;                                        % [deg] Leading edge radius
d = 0;                                              % [m] rocket tip-fin distance
zup_raw = 0.002;                                   % [m] fin semi-thickness
Lmaxu_raw = 0.002;                                  % [m] Fraction of chord from leading edge to max thickness


Xle1 = Lcenter + Lnose - d - C1;
diffC = C1-C2;

switch shape
    case 'rect'
        Xle2 = Lcenter + Lnose - d - C2;
        
    case 'iso'
        Xle2 = Lcenter + Lnose - d - diffC/2;
end

% Defining Fin Section
Zup = [zup_raw/C1, zup_raw/C2];
Lmaxu = [Lmaxu_raw/C1, Lmaxu_raw/C2];
Lflatu = [(C1 - 2*Lmaxu_raw)/C1, (C2 - 2*Lmaxu_raw)/C2];




%% Creating for005.dat file from previous data
if ismac   % mac procedure
    fid  =  fopen(strcat(pwd, '/for005.dat'),'w+');
else
    fid  =  fopen(strcat(pwd, '\for005.dat'),'w+');
end

%%%%%%%%%%%% Flight Conditions
%%%% Beta
fprintf(fid, '\n $FLTCON\r\n');
fprintf(fid, '  BETA = ');
fprintf(fid, '%.2f,\r\n', Beta(1));
%%%% Alt
fprintf(fid, ' ALT = ');
fprintf(fid, '%d', Nm);
fprintf(fid, '*');
fprintf(fid, '%d.,\r\n', Alt(1));
%%%% Nmach
fprintf(fid, '  NMACH = ');
fprintf(fid, '%d., \r\n', Nm);
%%%% Mach
fprintf(fid, '  MACH = ');
if Nm>11
    for M = 1:11
        fprintf(fid, '%.2f,',Mach(M));
    end
    fprintf(fid,  ' \r\n MACH(12) = ');
    for M = 12:Nm
        fprintf(fid, '%.2f',Mach(M));
        fprintf(fid, ',');
    end
else
    for M = 1:Nm
        fprintf(fid, '%.2f',Mach(M));
        fprintf(fid, ',');
    end
end
    
fprintf(fid, '\r\n');
%%%% Nalpha
fprintf(fid, '  NALPHA = ');
fprintf(fid, '%d., \r\n', Na);
%%%% Alpha
fprintf(fid, '  ALPHA = ');
for a = 1:9
    fprintf(fid, '%.1f,', Alpha(a));
end
fprintf(fid,  ' \r\n ALPHA(10) = ');
for a = 10:Na
    fprintf(fid, '%.1f,', Alpha(a));
end
fprintf(fid, '$');

%%%%%%%%%%%% Reference Quantities
fprintf(fid, '\r\n $REFQ\r\n');
%%%% xcg
fprintf(fid, '  XCG = ');
fprintf(fid, '%.4f,\r\n', xcg);
%%%% SREF
fprintf(fid, '  SREF = ');
fprintf(fid, '%.5f,\r\n', S);
%%%% LREF
fprintf(fid, '  LREF = ');
fprintf(fid, '%.3f,\r\n', D);
%%%% LATREF
fprintf(fid, '  LATREF = ');
fprintf(fid, '%.3f,$', D);

%%%%%%%%%%%% Axisymmetric Body Geometry
fprintf(fid, '\r\n $AXIBOD\r\n');
%%%% TNOSE
fprintf(fid, '  TNOSE = KARMAN, \r\n');
%%%% LNOSE
fprintf(fid, '  LNOSE = ');
fprintf(fid, '%.3f, \r\n', Lnose);
%%%% DNOSE
fprintf(fid, '  DNOSE = ');
fprintf(fid, '%.3f, \r\n', D);
%%%% LCENTR
fprintf(fid, '  LCENTR = ');
fprintf(fid, '%.3f, \r\n', Lcenter);
%%%% DCENTR
fprintf(fid, '  DCENTR = ');
fprintf(fid, '%.3f, \r\n', D);
%%%% DEXIT
fprintf(fid, '  DEXIT = ');
fprintf(fid, '%.3f, \r\n', D);
%%%% BASE
fprintf(fid, '  BASE = .FALSE.,$');

%%%%%%%%%%%% Finset
fprintf(fid, '\r\n $FINSET1 \r\n');
%%%% XLE
fprintf(fid, '  XLE = ');
fprintf(fid, '%.3f,', Xle1);
fprintf(fid, '%.3f, \r\n', Xle2);
%%%% NPANEL
fprintf(fid, '  NPANEL = ');
fprintf(fid, '%.1f,\r\n', Npanel);
%%%% PHIF
fprintf(fid, '  PHIF = ');
for P = 1:length(Phif)
    fprintf(fid, '%.1f', Phif(P));
    fprintf(fid, ',');
end
fprintf(fid, '\r\n');
%%%% LER
fprintf(fid, '  LER = 2*');
fprintf(fid, '%.4f,\r\n', Ler);
%%%% SSPAN
fprintf(fid, '  SSPAN = ');
fprintf(fid, '%.3f,', r);
fprintf(fid, '%.3f,\r\n', r + H);
%%%% CHORD
fprintf(fid, '  CHORD = ');
fprintf(fid, '%.3f,', C1);
fprintf(fid, '%.3f,\r\n', C2);
%%%% SECTYP
fprintf(fid, '  SECTYP = HEX, \r\n');
%%%% ZUPPER
fprintf(fid, '  ZUPPER = ');
fprintf(fid, '%.4f,', Zup(1));
fprintf(fid, '%.4f,\r\n', Zup(2));
%%%% LMAXU
fprintf(fid, '  LMAXU = ');
fprintf(fid, '%.4f,', Lmaxu(1));
fprintf(fid, '%.4f,\r\n', Lmaxu(2));
%%%% LMAXU
fprintf(fid, '  LFLATU = ');
fprintf(fid, '%.4f,', Lflatu(1));
fprintf(fid, '%.4f,$ \r\n', Lflatu(2));

%%%%%%%%%%%% Options
fprintf(fid, 'DERIV RAD \r\n');
fprintf(fid, 'DIM M \r\n');
fprintf(fid, 'DAMP \r\n');
fprintf(fid, 'SAVE \r\n');
fprintf(fid, 'NEXT CASE \r\n');

%%%%%%%%%%%% Cases
for A = 1:Nalt
    for B = 1:Nb
        if A == 1 && B == 1
        else
            fprintf(fid,' $FLTCON \r\n');
            fprintf(fid,'  BETA = ');
            fprintf(fid, '%.1f, \r\n', Beta(B));
            fprintf(fid, ' ALT = ');
            fprintf(fid, '%d', Nm);
            fprintf(fid, '*');
            fprintf(fid, '%d.,$ \r\n', Alt(A));
            fprintf(fid, 'DERIV RAD\r\n');
            fprintf(fid, 'DIM M\r\n');
            fprintf(fid, 'DAMP\r\n');
            fprintf(fid, 'SAVE\r\n');
            fprintf(fid, 'NEXT CASE\r\n');
        end
    end
end
fclose(fid);

%% Datcom and parsing
if ismac
    system('./datcom for005.dat' );
else
    system('datcom.exe for005.dat' );
end


value = 0;
while value == 0
    value = exist('for006.dat','file');
    pause(0.1);
end
clc

[Coeffs, State] = datcom_parser(s.mat_name);

delete('for003.dat', 'for004.dat', 'for005.dat',  'for009.dat',...
    'for010.dat', 'for011.dat', 'for012.dat')


