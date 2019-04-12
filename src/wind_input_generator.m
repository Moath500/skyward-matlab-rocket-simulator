function [uw,vw,ww] = wind_input_generator(settings,z,uncert)

% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
% Release date: 13/01/2018

% This function allows to use a custom set of wind, defined in config.m 

magn = (1 + uncert(1)/100).*settings.wind.input_matr(1,:);
dir = mod(180 + settings.wind.input_matr(2,:),360);
dir = dir + uncert(2);

uw_vect = magn.*cosd(dir);
vw_vect = magn.*sind(dir);
h_vect = settings.wind.input_matr(3,:);

h = -z;

if h < 0
    h = 0;
end

uw = interp1(h_vect,uw_vect,h);
vw = interp1(h_vect,vw_vect,h);
ww = 0;