clear
close all
clc;
%%
dataPath = '../dataset_annotated'; 
visPath = '../visAnnotation';
destPath = '../dataSet_wormBody'; % save the result under this folder
dataPathStorage = '.\dataset_segmentation';
%%
imList = dir(fullfile(destPath, '*.mat'));
for j = 1:length(imList)
    load(fullfile(destPath, imList(j).name));
    
    for i = 1:length(wormSetMore)
        wormSetMore{i}.im = [];
    end
    
    save(fullfile(destPath, imList(j).name), 'mask', 'wormSetMore' );
end
%% leaving blank


