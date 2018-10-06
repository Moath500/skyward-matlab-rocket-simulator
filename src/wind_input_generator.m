function [uw,vw,ww] = wind_input_generator(settings,z)
% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
% Release date: 13/01/2018

% This function allows to use a custom set of wind, defined in config.m 

uw_vect = settings.wind.input_matr(1,:).*cosd(360-settings.wind.input_matr(2,:));
vw_vect = settings.wind.input_matr(1,:).*sind(360-settings.wind.input_matr(2,:));
h_vect = settings.wind.input_matr(3,:);

h = -z+settings.z0;
if h < 0
    h = 0;
end

uw = interp1(h_vect,uw_vect,h);
vw = interp1(h_vect,vw_vect,h);
ww = 0;