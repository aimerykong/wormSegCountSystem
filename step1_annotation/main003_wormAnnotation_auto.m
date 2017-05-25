clear
clc
close all;
%% configuration
curImgId = 7;
flag_using_histeq = true; % true false
numWorm = 10;  % the number of worms required to find by users

radius = 40;

%% annotating the worms
extractWormBody_interface_V5(curImgId, flag_using_histeq, numWorm, radius );