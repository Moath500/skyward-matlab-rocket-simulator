function [ww] = vert_windgen(MagMin,MagMax)
% vert_windgen(MagMin,MagMax)
%function that generates vertical wind component

% Author: Ruben Di Battista & Gabriele Poiana
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu gabriele.poiana@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 17.I.2016
% License:  2-clause BSD

%Generating random value for magnitude

x = rand;
ww = MagMin+(MagMax-MagMin)*x;
if x < 0.5
    ww = -ww; 
end
end

