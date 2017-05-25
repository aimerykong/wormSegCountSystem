clear all force
clc
close all;
%% configuration
curImgId = 97;
numWorm = 5;  % the number of worms required to find by users

flag_using_histeq = false; % true false

radius = 80;
flag_using_manual = false; % true false
dataPath = '../dataset_annotated';
%% annotating the worms
if flag_using_manual % manually drawing the worm body
    extractWormBody_interface_manual(curImgId, flag_using_histeq, numWorm, radius, dataPath );
else % semi-automatic annotation using dynamic chain model, allowing for correcting worm manually
    extractWormBody_interface_auto(curImgId, flag_using_histeq, numWorm, radius );    
end
%% leaving blank

