clear all
close all
clc

D = 0.09; %Diameter [m]
L = 2.02; %Length [m]

Chord_1 = 0.14; %[m]
Chord_2 = 0.07; %[m]
Alt = 0.07; %[m]

Spess = 0.004; %[m]
Dist = 0.00; %[m]
Ler = 0.002;


XLE_1 = L-Chord_1-Dist;
XLE_2 = L-Chord_2-Dist; %Trapezioidali

%XLE_2 = L-Chord_1+(Chord_1-Chord_2)/2-Dist;

SSPAN_1 = D/2;
SSPAN_2 = D/2 + Alt;

ZUPPER_1 = (Spess/2)/Chord_1;
ZUPPER_2 = (Spess/2)/Chord_2;

LMAXU_1 = Ler/Chord_1;
LMAXU_2 = Ler/Chord_2;

LFLATU_1 = (Chord_1 - 2*Ler)/Chord_1;
LFLATU_2 = (Chord_2 - 2*Ler)/Chord_2;

fprintf('XLE = %f , %f\n\n',XLE_1,XLE_2);
fprintf('SSPAN = %f , %f\n\n',SSPAN_1,SSPAN_2);
fprintf('ZUPPER = %f , %f\n\n',ZUPPER_1,ZUPPER_2);
fprintf('LMAXU = %f , %f\n\n',LMAXU_1,LMAXU_2);
fprintf('LFLATU = %f , %f\n\n',LFLATU_1,LFLATU_2);
