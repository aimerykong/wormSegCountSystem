function [wormSetMore] = rotflip2getMore(imOrg, wormFound)
% generate more annotated image by rotating the original image with its
% worm labels.
%
%
% see also compact2full, showFoundWorm
%
% Shu Kong
% 04/13/2015

%% generate more annotated image by rotating the original image with its worm labels
wormSetMore = cell(8,1);

count = 1; % original
wormSetMore{count}.im = imOrg;
wormSetMore{count}.wormFound = wormFound;


count = count + 1; % rotate counterclockwise 90 degree
wormSetMore{count}.im = 'rot90(imOrg, 1)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = size(imOrg,1) - wormFound{i}(2,:);
    wormFound2{i}(2,:) = wormFound{i}(1,:); 
end
wormSetMore{count}.wormFound = wormFound2;
        

count = count + 1; % rotate counterclockwise 180 degree
wormSetMore{count}.im = 'rot90(imOrg, 2)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = size(imOrg,2) - wormFound{i}(1,:);
    wormFound2{i}(2,:) = size(imOrg,1) - wormFound{i}(2,:);
end
wormSetMore{count}.wormFound = wormFound2;
    

count = count + 1; % rotate counterclockwise 270 degree
wormSetMore{count}.im = 'rot90(imOrg, 3)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = wormFound{i}(2,:);
    wormFound2{i}(2,:) = size(imOrg,2) - wormFound{i}(1,:);
end
wormSetMore{count}.wormFound = wormFound2;


%% flip left right and do the same thing
count = count + 1; % flip left right
imOrg = fliplr(imOrg);
for i = 1:numel(wormFound)
    wormFound{i}(2,:) = size(imOrg,2) - wormFound{i}(2,:);
end
wormSetMore{count}.im = 'fliplr(imOrg)';%imOrg;
wormSetMore{count}.wormFound = wormFound;
    

count = count + 1; % rotate LR-flipped image counterclockwise 90 degree
wormSetMore{count}.im = 'rot90(imOrg, 1)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = size(imOrg,1) - wormFound{i}(2,:);
    wormFound2{i}(2,:) = wormFound{i}(1,:); 
end
wormSetMore{count}.wormFound = wormFound2;
        

count = count + 1; % rotate LR-flipped image counterclockwise 180 degree
wormSetMore{count}.im = 'rot90(imOrg, 2)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = size(imOrg,2) - wormFound{i}(1,:);
    wormFound2{i}(2,:) = size(imOrg,1) - wormFound{i}(2,:);
end
wormSetMore{count}.wormFound = wormFound2;
    

count = count + 1; % rotate LR-flipped image ounterclockwise 270 degree
wormSetMore{count}.im = 'rot90(imOrg, 3)';
wormFound2 = wormFound;
for i = 1:numel(wormFound)
    wormFound2{i}(1,:) = wormFound{i}(2,:);
    wormFound2{i}(2,:) = size(imOrg,2) - wormFound{i}(1,:);
end
wormSetMore{count}.wormFound = wormFound2;

    









