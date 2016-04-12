close all
clear, clc

[~, a, ~, rho] = atmoscoesa(0);

% DATCOM DATA
load for006_full.mat;
CoeffsF = Coeffs;
load for006_empty.mat;
CoeffsE = Coeffs;

givA = State.Alphas;
givB = State.Betas;
givH = State.Altitudes;
givM = State.Machs;

% GEOMETRIC DATA
C = 0.174; %m Reference Length (Fuselage Diameter)
S = 0.02378; %m2 Cross-sectional Surface
L = 4.30; %m lunghezza missile

% DRAG CURVE
M = linspace(0, 3);
V = a*M;
CD = zeros(1, length(M));
D = zeros(1, length(M));
for n = 1:length(M)
    CD(n) = interpn(givA,givM,givB,givH,CoeffsF.CD,0,M(n),0,0);
    D(n) = 1/2 * rho * V(n)^2 * S * CD(n);
end

figure();
subplot(2,1,1)
plot(V, CD), grid on;
ylabel('CD [-]')
xlabel('V [m/s]')
ylim([0, 0.7])

% figure();
subplot(2,1,2)
plot(V, D), grid on;
ylabel('D [N]')
xlabel('V [m/s]')