% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

clear, clc

% %Cesaroni PRO 150 White Thunder
% sTH = [0 0.05 0.15 0.5 0.6 0.74 0.85 1.15 1.35 2.45 ...
%     3 3.7 4 4.5 4.8 4.9 5 5.05 5.1 5.15 5.2]; %s
% THs = [0 9200 7900 8500 8500 8350 8300 8400 8500 8500 ...
%     8300 8000 7800 7600 7450 7300 7300 4500 500 100 0]; %N

%Cesaroni PRO 150 SkidMark
sTH = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.2 1.8 3.2 3.6 ...
    4.8 6 7 7.2 7.6 7.8 7.9 8 8.1 8.19]; %s
THs = [0 3400 3100 3000 3300 3400 3500 3700 3700 3800 ...
    4000 4081.6 3900 3800 3700 3500 3350 3200 3000 2000 750 0]; %N

s=linspace(0,sTH(end));

dt = s(2);
TH = interp1(sTH,THs,s);

sum = 0;
for j = 1:length(s)
    sum = sum + TH(j)*dt;
end
media = sum/5.12 %N

plot(s,TH);
grid on;
