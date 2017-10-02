clear; clc; close all

%% EARTH
% no rogallo
figure();
back('circumference.fig',[153/255 153/255 0],[1 1 204/255])

% rogallo
figure();
back('triangle.fig',[0 102/255 0],[204/255 1 204/255])

%% SEA
%rogallo
figure();
front('triangle.fig',[0 102/255 204/255],[204/255 1 1])

%ballistic (45 deg)
figure();
cone('triangle.fig', [1 229/255 204/255])

%apogee
figure();
cone('triangle.fig', [224/255 224/255 224/255])