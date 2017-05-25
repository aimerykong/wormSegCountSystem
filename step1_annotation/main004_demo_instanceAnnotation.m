% ---- worm dataset label generation (version 1.0) ----
% Given the worm body lines, this code generates segmentation map, keypoint
% heatmap, and bounding boxes for detection.
%
%
% Through analyzing the original images, we find the longest worm occupies 
% ~220 pixels, while the smallest one takes about 40 pixel. We downsize the
% image to be of half size, then worms are in the range of [20, 110]
% pixels for the longer side, either width or length.
%
% We are interested in borrowing the VGG16 model, where at layer conv5_3 in
% VGG16 the receptive field size reflected in the input image space is
% 196x196. Therefore, we choose to use layers before, including, conv5_3
% VGG16. There are four max-pooling layer, therefore the heatmap is 4 times
% smaller than the input image.
%
% We fine-tune the model and train additional layers on top of
% conv5_3 to predict the "label", e.g. segmentation map, keypoint heatmap. 
%
% Segmentation map can be used as regularizer.
%
% The "label" should be vector for each pixel, including
%   [conf of keypoint;
%    width-offset; 
%    height-offset ]
%
%
% Shu Kong
% skong2@ics.uci.edu
% 05/22/2017
clear
close all
clc;

dataPath = '../dataset_annotated'; 
visPath = '../visAnnotation';
destPath = '../dataSet_wormBody'; % save the result under this folder
dataPathStorage = '.\dataset_segmentation';

imList = dir(fullfile(dataPath, '*.png'));

%% parameters and path
curImgId = 40;

filename = imList(curImgId).name; fprintf('\n%s\n',filename);
[junk,nameFile,extFile] = fileparts(filename); % get the file name
%% imread and display for clicking
im = imread( fullfile(dataPath, filename) ); % read image
% im = im(1:round(size(im,1)*0.7), round(size(im,1)*0.1):round(size(im,1)*0.8)); % crop the image to a reasonable area
imOrg = im; % backup the original image

dirMat = dir( fullfile(destPath, strcat(nameFile, '*.mat')) ); % retrieval all the existing images

wormListAll = {};
mask = single(0*im);

[orgSizeY, orgSizeX] = size(mask); % mask of original image
for i = 1:numel(dirMat) % merge all the existing masks to visualize how many worms are still needed to estimate
    matTMP = load( fullfile(destPath, dirMat(i).name) ); % load the mask stored in the directory
    for j = 1:numel(matTMP.wormSetMore{1}.wormFound)
        wormListAll{end+1} = matTMP.wormSetMore{1}.wormFound{j};        
    end
end
%% ---------- show image ----------------------------------------
figure(1);
subplot(1,2,1);
imshow(imOrg); title('original image'); 
for i = 1:numel(wormListAll) % draw worm body on the original image
    A = xor(mask, single(connPart4Body(mask, wormListAll{i})));
    mask = mask + (i*single(A));
    line( wormListAll{i}(2,1:end), wormListAll{i}(1,1:end), 'linewidth', 1, 'color', 'r');    
end
%% ---------- show segmentation mask ----------------------------------------
subplot(1,2,2);
SE = strel('disk', 2);
mask = imdilate(single(mask), SE);
imagesc(mask); axis image off; %  colorbar;
title('label mask');
%% leaving blank

