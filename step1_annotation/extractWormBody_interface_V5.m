function extractWormBody_interface_V5(curImgId, flag_using_histeq, numWorm, radius)
% ---- worm dataset collection (version 1.0) ----
% click the worm head and tail (body parts) to estimate worm body
% -------------------------------------------------
% With the worm head and tail, the body can be estimated by a chain
% model with dynamic programming. While when this method fails, the user
% can choose to draw a line/curve to annotate the worm.
%
% There are multiple configuration choices for a user to set --
%   1. the path of the folder which contains the images to annotate;
%   2. the name of the current image to annotate;
%   3. whether to use histogram equalization to augment the intensity
%       contrast;
%   4. the number of worms to annotate for the current run;
%   5. others including the foler to store annotation, visualization, etc.
% 
%
% To annotate the worm, first click a body, preferrably at the middle of
% worm body, then a zoom-in image will be displayed for clicking the worm
% head and tail. After that, the estimated worm will be displayed to
% indicate how well the method estimates the worm.
% If the result is good enough, the user can submit the result. If the
% result is not good enough, the user can choose either to return to the
% original image and start to annotate a new worm, or manually draw a curve
% to annotate the worm.
%
% An advice is to try this UI with personal expectations before annotating
% formally.
%
% Thanks to Sahand Ghorbanpour for his effort in making the UI more
% approachable.
%
%
% Shu Kong
% skong2@uci.edu
% 05/16/2017

%% parameters and path
% curImgId = 2;
% flag_using_histeq = false; % true false
% numWorm = 2;  % the number of worms required to find by users


dataPath = '../dataset_annotated'; % the images to estimate worms are under this folder
imList = dir(fullfile(dataPath, '*.png'));

% filename = 'scan2_plate9.png';% 
filename = imList(curImgId).name; fprintf('\n%s\n',filename);

% radius = 80; % the patch size is of (2*radius+1)x(2*radius+1) pixel resolution
scoreMultiplier = 20; % calibrating the mismatching score
nPart = 12; % the number of parts for chain model to estimate worm body
numParts = nPart+2; % click 10 parts for one worm
partSize = 3; % define the part size % 5
subSize = 2; % define the subsampling straddle % 3
partDist = 3; % partwise distance % 4

% radius = 80; % the patch size is of (2*radius+1)x(2*radius+1) pixel resolution
% nPart = 14; % the number of parts for chain model to estimate worm body
% numParts = nPart+2; % click 10 parts for one worm
% partSize = 5; % define the part size % 5
% subSize = 3; % define the subsampling straddle % 3
% partDist = 6; % partwise distance % 4
% numPartMulipleClicking = 10; % the number of clicks required to estimate the worm manually (clicks at worm parts in order)

visPath = '../visAnnotation';
destPath = '../dataSet_wormBody'; % save the result under this folder
%% imread and display for clicking
im = imread( fullfile(dataPath, filename) ); % read image

if ~isdir(visPath)
    mkdir(visPath);
end

% im = im(1:round(size(im,1)*0.7), round(size(im,1)*0.1):round(size(im,1)*0.8)); % crop the image to a reasonable area
imOrg = im; % backup the original image
if flag_using_histeq
    im = histeq(im);
end

mask = zeros(size(imOrg)); % the mask as labels for later use, will be saved along with the estimated worms

imDisplay = repmat(im, [1,1,3]); % to display the estimated worms in the whole original image
maskDisplay = mask; % to display the estimated worms in binary scale with original size

[~, nameFile, ~] = fileparts(filename); % get the file name
dirMat = dir( fullfile(destPath, strcat(nameFile, '*.mat')) ); % retrieval all the existing images

wormListAll = {};
for i = 1:numel(dirMat) % merge all the existing masks to visualize how many worms are still needed to estimate
    matTMP = load( fullfile(destPath, dirMat(i).name) ); % load the mask stored in the directory
    maskDisplay = maskDisplay | matTMP.mask; % 'or' operation
    for j = 1:numel(matTMP.wormSetMore{1}.wormFound)
        wormListAll{end+1} = matTMP.wormSetMore{1}.wormFound{j}; %#ok<AGROW>
    end
end

% display the image in rgb format for better visualization
imDisplay = getOverallMask(imDisplay, maskDisplay);

imOrg = imDisplay;

%% procedure of estimating the worms
wormFound = cell(numWorm,1); % store all the found worms by the user
posWorm = zeros(numWorm, 7); % one row per worm, [center_x, center_y, upperleft_x, upperleft_y, bottomright_x, bottomright_x, num]

countWorm = 0;
retakeFlag = 0; % flag to re-estimate the worm
tenPartClickFlag = 0; % flag to indicate whether the user wants to click ten times manually to estimate the worms
manualDrawDoneFlag = 0;


while countWorm < numWorm % get all valid worms
    if ~retakeFlag % if user does not want to re-estimate the previous worm, then click the worm to zoom in for estimation
        %figure(1); showFoundWorm(imOrg, wormFound, countWorm);
        figure(1);
        showFoundWorm(imOrg, wormListAll, numel(wormListAll));
        %delete(findall(findall(gcf,'Type','axe'),'Type','text'))
        title(strcat(num2str(numWorm-countWorm), ' worms are waiting for you by clicking their head and tail'))
        wormPartPosition = zeros(2, numParts);
        
        [xPatchCenter,yPatchCenter] = ginput(1); % zoom in the area centered at (x, y)
    end
    
    %make some arbitrary rectangle (in this case, located at (0,0) with [width, height] of [10, 20])
    UL = [yPatchCenter, xPatchCenter]-radius;
    BR = [yPatchCenter, xPatchCenter]+radius;
    UL(UL<=0) = 1;
    
    if BR(1) > size(imOrg, 1) % please make sure all the images are square patch, do not click those which are too close to the image boundary
        BR(1) = size(imOrg, 1);     end
    if BR(2) > size(imOrg, 2)
        BR(2) = size(imOrg, 2);    end
    ULorg = int16(UL); % convert to int for index reason
    BRorg = int16(BR);
    
    patchTMP = imOrg(ULorg(1):BRorg(1), ULorg(2):BRorg(2)); % get the zoom-in area
    
    patchTMP = repmat(patchTMP, [1,1,3]); % make it a rgb format for annotation reason
    S.fH = figure(1);
    %     imshow(patchTMP); title('click worm head and tail');
    S.aH = axes;
    S.iH = imagesc(patchTMP); 
    delete(findall(findall(gcf,'Type','axe'),'Type','text'))
    title('click worm head and tail');
    patchTMP_BACKUP = patchTMP;
    
    
    if tenPartClickFlag % if the user chooses to click ten times to estimate the worm manually
        %{
        countParts = 1;
        wormPartPosition = zeros(2,numPartMulipleClicking);
        while countParts <= numPartMulipleClicking % click the worm parts and estimate the worm body based on the parts
            [x,y] = ginput(1);
            wormPartPosition(:, countParts) = [y;x];
            
            patchTMP(uint8(y), uint8(x), 1:3) = [255, 0, 0]; % set the clicked part points as red
            figure(1); imshow(patchTMP); title( strcat(num2str(countParts), '/', num2str(numPartMulipleClicking), ' click worm body in order') );
            line( wormPartPosition(2,1:countParts), wormPartPosition(1,1:countParts), 'linewidth', 1, 'color', 'r');
            
            countParts = countParts + 1;
        end
        %}
        hold on
        X = [];
        Y = [];
        set(S.iH,'ButtonDownFcn',@startDragFcn)
        set(S.fH, 'WindowButtonUpFcn', @stopDragFcn);
        delay = 0.01;  % 10 milliseconds
        while ~manualDrawDoneFlag  % set by the callback
            pause(delay);  % a slight pause to let all the data gather
        end
        
        
        wormPartPosition(1,:) = double(ULorg(1)) + wormPartPosition(1,:); % y -- vertical (1)
        wormPartPosition(2,:) = double(ULorg(2)) + wormPartPosition(2,:); % x -- horizontal (2)
        subsampleX = 1:3:size(X,2);
        subsampleY = 1:3:size(Y,2);
        X = X(:,subsampleX);
        Y = Y(:, subsampleY);
        wormPartPosition = [Y;X];
        manualDrawDoneFlag = 0;
        % line( X, Y, 'linewidth', 2,'LineStyle',':');
        hold off;
        tenPartClickFlag = 0;
        
    else % the user click only the head and tail of the worm for automatically estimating the worm
        countParts = 1;
        while countParts <= 2 % click the worm parts and estimate the worm body based on the parts
            [x,y] = ginput(1);
            disp(x);
            disp(y);
            wormPartPosition(:, countParts) = [y;x];
            patchTMP(uint16(y), uint16(x), 1:3) = [255, 0, 0]; % set the clicked part points as red
            countParts = countParts + 1;
        end
        
        % dynamic programming to estimate the worm body
        wormPartPositionBACKUP = wormPartPosition;
        %  make clicked head and tail set to 1 and all other bits to zero
        target = zeros(size(patchTMP,1), size(patchTMP,2) );% which dimension = 1 or 2
        target( sub2ind(size(target), uint16(wormPartPositionBACKUP(1,1)),uint16(wormPartPositionBACKUP(2,1))) ) = 1;
        target( sub2ind(size(target), uint16(wormPartPositionBACKUP(1,2)),uint16(wormPartPositionBACKUP(2,2))) ) = 1;
        
        % threshold the patch to faciliate automatic estimation
        p = 20; % p-th percentile of the intensities in the patch as the threshold
        thresh = patchTMP(:); % flatten the patch Array
        thresh = prctile(thresh, p); % get the threshold intensity value
        patchTMP_BACKUP(patchTMP_BACKUP >= thresh) = Inf; % set all pixels greater than threshod to inf
        
        % chain model via dynamic programming to estimate the worm body
        wormStruct = dpWormEstimation(patchTMP, mat2gray(patchTMP_BACKUP(:,:,1)), target, partDist, scoreMultiplier, partSize, subSize, nPart );
        % patchTMP(uint16(wormStruct.parts(:,2)), uint16(wormStruct.parts(:,1)),1:3) = [255, 0, 0];
        for i = 1:wormStruct.size
            patchTMP(wormStruct.parts(i,1),wormStruct.parts(i,2), 1:3) = [255, 0, 0];
        end
        %  patchTMP(wormStruct.parts(3,1),wormStruct.parts(3,2), 1:3) = [255, 0, 0];
        %  patchTMP(wormStruct.parts(4,1),wormStruct.parts(4,2), 1:3) = [255, 0, 0];
        % patchTMP(wormStruct.parts(5,1),wormStruct.parts(5,2), 1:3) = [255, 0, 0];
        %patchTMP(wormStruct.cx,wormStruct.cy, 1:3) = [255, 0, 0];
        % end
        figure(1);
        %         imshow(patchTMP); title('estimated worm body');
        imagesc(patchTMP); 
        delete(findall(findall(gcf,'Type','axe'),'Type','text'))
        title('estimated worm body');
        % set the clicked part points as red
        line( wormStruct.parts(:,2), wormStruct.parts(:,1), 'linewidth', 3,'LineStyle',':');
        
        
        wormPartPosition = wormStruct.parts(:,1:2)';
    end
    wormPartPosition(1,:) = double(ULorg(1)) + wormPartPosition(1,:); % y -- vertical (1)
    wormPartPosition(2,:) = double(ULorg(2)) + wormPartPosition(2,:); % x -- horizontal (2)
    
    % get the mask as labels
    mask = connPart4Body(mask, wormPartPosition);
    
    % maskDisplay = maskDisplay | mask; % store into the overall mask
    imDisplay = getOverallMask(imDisplay, maskDisplay);
    
    ButtonName = questdlg('Are you sure to submit? Please make sure you estimate the worm body correctly! \nYes: submit, No: do not submit and turn to another patch, ... Re-take: re-click the worm on this patch', ...
        'double check before submission', ...
        'Yes', 'No (return)', 'Manual Drawing',  'No');
    switch ButtonName,
        case 'Yes', % if yes, store all the patches
            retakeFlag =0;
            tenPartClickFlag = 0;
            countWorm = countWorm + 1;
            wormFound{countWorm} = wormPartPosition;
            wormListAll{end+1} = wormPartPosition;
            
            maskDisplay = maskDisplay | mask; % store into the overall mask
            imDisplay = getOverallMask(imDisplay, maskDisplay);
            
            %showFoundWorm(imOrg, wormFound, countWorm);
            disp('Thank you for your contribution, please do more!');
        case 'Manual Drawing',
            retakeFlag = 1;
            tenPartClickFlag = 1;
            delete(findall(findall(gcf,'Type','axe'),'Type','text'))
        case 'No (return)',
            retakeFlag = 0;
            tenPartClickFlag = 0;
            disp('Thank you for your carefulness, please redo it!');
    end
    clf;
end

%% flip and rotate to get more annotated data
[wormSetMore] = rotflip2getMore(imOrg, wormFound);

if ~isdir(destPath)
    mkdir(destPath);
end

%% save the result
[~, nameFile, ~] = fileparts(filename);
save( fullfile(destPath, strcat(nameFile, '_', datestr(now,'mm-dd-yyyy-HH-MM-SS'), '__', num2str(numWorm), '.mat')),...
    'wormSetMore', 'mask');%
 
%% final showcase the annotated worms
wormSetMore = compact2full(wormSetMore);
%k = 1; showFoundWorm(wormSetMore{k}.im, wormSetMore{k}.wormFound);

f = figure(1);
showFoundWorm(im, wormListAll, numel(wormListAll));
H = getframe(f);
delete(findall(findall(gcf,'Type','axe'),'Type','text'))
title('visualization of estimated worms');

[~,newfilename,~] = fileparts(filename);
print(f, '-r100', '-dbitmap', fullfile(visPath, [newfilename '_vis.bmp']));


    function startDragFcn(varargin)
        set( S.fH, 'WindowButtonMotionFcn', @draggingFcn );
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);
        X = x;
        Y = y;
    end


    function draggingFcn(varargin)
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);
        X = [X x];
        Y = [Y y];
        plot(X,Y,'r','LineWidth',4,'ButtonDownFcn',@startDragFcn)
        hold on
        drawnow
    end

    function stopDragFcn(varargin)
        set(S.fH, 'WindowButtonMotionFcn', '');  %eliminate fcn on releas
        if tenPartClickFlag
            manualDrawDoneFlag = 1;
            % line( X, Y, 'linewidth', 10,'LineStyle',':');
        end
    end
end

