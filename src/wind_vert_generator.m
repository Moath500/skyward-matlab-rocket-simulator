function [ww] = wind_vert_generator(MagMin, MagMax)
%{

wind_vert_generator - function that generates vertical wind component

INPUTS:
            - MagMin, Minimum wind magnitude;
            - MagMax, Maximum wind magnitude.

OUTPUTS:
            - ww, wind component along z.

Author: Ruben Di Battista & Gabriele Poiana
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu gabriele.poiana@skywarder.eu
April 2014; Last revision: 17.I.2016

%}

%Generating random value for magnitude

x = rand;
ww = MagMin + (MagMax - MagMin)*x;
if x < 0.5
    ww = -ww; 
end
end

