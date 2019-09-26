close all
clear all
clc
direction = 192;
ramp_direction = 200;
real_direction = ramp_direction + direction -180 ;

if real_direction > 360
    direction = real_direction - 360;
end

speed = 3.2;
speed = speed/3.6
direction