clear
clc
close all;
%% configuration
curImgId = 20;
numWorm = 2;  % the number of worms required to find by users

flag_using_histeq = false; % true false

radius = 70;

%% annotating the worms
extractWormBody_interface_manual(curImgId, flag_using_histeq, numWorm, radius );