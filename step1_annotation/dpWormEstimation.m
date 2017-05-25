function [wormStruct] = dpWormEstimation(imOrg, im, target, partDist, scoreMultiplier, partSize, subSize, nPart )
% This function is used for estimating worm body with a simple part-based
% model, which is a chain model. Dynamic programming is used to find a
% global optimal placement of the worm parts. This function returns the
% worm parts locations and estimate score stored in a structure.
%
%
% Input
%   im              --  the image or image area;
%   target          --  the (predicted) target map highlighting the TWO 
%                       detected keypoints of the worms. It is a binary
%                       map.
%   partDist        --  partwise distance as a defined parameter, e.g. 3 as
%                       default;
%   scoreMultiplier --  calibrating the mismatching score into the same
%                       scale of partwise distance, e.g. 10 as default;
%   partSize        --  define the part size, e.g. 5x5 as default; 
%   subSize         --  define the subsampling straddle, e.g. 3 as default;
%
%
% Output
%   wormStruct      --  a structure with part locations, total mismatch
%                       score (the smaller the better).
%
% Shu Kong @ UCI
% skong2@uci.edu
% 03/21/2015

%% setting default values
if (nargin < 3) 
    partDist = 3; end
if (nargin < 4) 
    scoreMultiplier = 10; end
if (nargin < 5) 
    partSize = 5; end
if (nargin < 6) 
    subSize = 5; end
if (nargin < 7) 
    nPart = 14; end % number of parts

%% build fast intermediate maps for dynamic programming
keypointIdx = find(target==1); % find the two detected keypoints
[keypointsubI, keypointsubJ] = ind2sub(size(im), keypointIdx); % subscripts

% define average filter for fast mismatching calculation
h = fspecial('average', partSize); 
imMeanFF = filter2(h, im, 'same');
imMeanFF(1:ceil(partSize/2),:) = 1; % boundary does not contribute at all
imMeanFF(:, 1:ceil(partSize/2)) = 1;
imMeanFF(end-ceil(partSize/2):end,:) = 1;
imMeanFF(:, end-ceil(partSize/2):end) = 1;

% subsample the image into a smaller one, for fast calculation and saving memory
subsampleRow = 1:subSize:size(im,1); 
subsampleColumn = 1:subSize:size(im,2);
imMeanF = imMeanFF(subsampleRow, subsampleColumn);

% get the subscript and linear index in the subsampled image
keypointsubI = ceil(keypointsubI/subSize);
keypointsubJ = ceil(keypointsubJ/subSize);
keypointIdxNew = sub2ind(size(imMeanF), keypointsubI, keypointsubJ);

C = imMeanF*scoreMultiplier; % store the mismatch score for fast calculation
V = zeros(size(C,1), size(C,2), numel(C)); % storing the pairwise distances
for i = 1:size(V,3) % fast index along the third mode for partwise distance
    D = reshape(1:numel(C), [], 1);
    [I, J] = ind2sub(size(C), D);
    [anchorI, anchorJ] = ind2sub(size(C), i);
    I = I-anchorI;
    J = J-anchorJ;
    D = (I.^2+J.^2);
    
    %D = max(0, (D-partDist).^2-4).^3;
    %D = (sqrt(I.^2+J.^2)-partDist).^2;
    %D = abs(sqrt(I.^2+J.^2)-partDist).^2; % distance can be chosen to be other forms, say Manhattan distance
                                          % 'partDist' is the desired
                                          % distance.
    V(:,:,i) = reshape(D, size(C));
end

%% dynamic programming
L0 = keypointIdxNew(1); % worm head
Ln1 = keypointIdxNew(2); % worm tail

B = zeros(nPart, size(V,3)); % the table used for dynamic programming
V = reshape(V, numel(V(:,:,1)), size(V,3));
for i = 1:nPart
    if i == 1
        TMP = V( L0, : );
        B(i,:) = (C(:) + TMP(:))';
    else
        %B(i,li) = C(li) + min(TMP) + (i==nPart)*V(li, Ln1); % the recurrent calculation
        TMP = bsxfun(@plus, V, B(i-1,:) );
        B(i,:) = C(:) + min(TMP, [], 2) + (i==nPart)*V(:, Ln1);
    end
end

%% trace back to find the part sequence
parts = zeros(1, nPart+2); parts(1) = Ln1;
mask = zeros(size(C));
mask(Ln1) = 1;

imDemons = imMeanF; % for demonstration reason
imDemons(keypointIdxNew) = 1.2;
[score, ind] = min(B(nPart, :));
imDemons(ind) = 1.2;
mask(ind) = 1;
parts(2) = ind; % the last part excluding the worm tail
[ex,ey] = ind2sub(size(imMeanF),Ln1);
[sx,sy] = ind2sub(size(imMeanF),L0);
[xs,ys] = ind2sub(size(imMeanF),ind);
val=imMeanF(xs,ys);

%disp(val);

for i = nPart:-1:1+1 
    li = ind;
    [val, ind] = min( V(li,:) + B(i-1,:) );
    imDemons(ind) = 1.2;
    mask(ind) = 1;
    parts(nPart-i+3) = ind;
%     fprintf('\t%d\n', ind);
end
mask(L0) = 1; % worm head
parts(end) = L0;
threshold = 0.65;
if  val >  0.1
   %point =  findNeighbor(imMeanF,xt,yt,3,threshold);
   %val=imMeanF(point.x,point.y);
   %parts(2) = sub2ind(size(imMeanF), point.x, point.y);
   result = findConnectionParts(imMeanF,sx,sy,ex,ey,1,threshold);
   %[xParts,yParts] = ind2sub(size(imMeanF), result.points);
   xParts = result.points(:,1);
   yParts = result.points(:,2);
   xParts = (xParts-1)*subSize+1;
   yParts = (yParts-1)*subSize+1;
end 
    

%% get the parts location with the original image size
%[xParts,yParts] = ind2sub(size(mask), parts);
%xParts = (xParts-1)*subSize+1;
%yParts = (yParts-1)*subSize+1;

mask = zeros(size(im,1), size(im,2)); %imresize(0*mask, size(mask)*subSize);
idx = sub2ind(size(mask), xParts, yParts);
mask(idx) = 1;

%% assemble the output
wormStruct.size=size(result.points,1);
wormStruct.showMap = imDemons;
wormStruct.mask = mask;
wormStruct.parts = [xParts(:), yParts(:)];
wormStruct.score = score + ...
    1-target(L0) + ...
    1-target(Ln1);



function [result]= findConnectionParts(im,StartX,StartY,EndX,EndY,maskDim,threshold)
    converged = 0;
    maxIteration = 40;
    pointS = zeros(1,maxIteration);
    pointE = zeros(1,maxIteration);
    sfoundPoints = zeros(1,maxIteration);
    efoundPoints = zeros(1,maxIteration);
    oldSX = intmax;
    oldSY = intmax;
    oldEX = intmax;
    oldEY = intmax;
    headOrTail=true;
    pointS(1,1)= sub2ind(size(im), StartX, StartY);
    pointE(1,1)= sub2ind(size(im), EndX, EndY);
    sfoundPoints(1,1) =  pointS(1,1);
    efoundPoints(1,1) =  pointE(1,1);
    sIt = 2;  % iterator
    eIt = 2;
    it =  1;
    


    mapPoints = containers.Map('KeyType','int32','ValueType','int32');
    while  it <= maxIteration % to get sure not trapped  in infinity if they don't converge
        nextNeighborH = findNeighbor(im,StartX,StartY,maskDim,threshold,oldSX,oldSY,headOrTail,sfoundPoints,mapPoints); % starting from Head
        nextNeighbotT = findNeighbor(im,EndX,EndY,maskDim,threshold,oldEX,oldEY,headOrTail,efoundPoints,mapPoints); % staring  from Tail
        intensityS    =  nextNeighborH.intensity;
        intensityE    =  nextNeighbotT.intensity;
     

        nextNeighborH = adjustInvalidPoint(nextNeighborH,maskDim);
        nextNeighbotT = adjustInvalidPoint(nextNeighbotT,maskDim);  
        convertedIndexS = sub2ind(size(im), nextNeighborH.x, nextNeighborH.y);
        convertedIndexE =  sub2ind(size(im), nextNeighbotT.x, nextNeighbotT.y);
        
        if intensityS > threshold
            mapPoints(convertedIndexS) = intensityS;
        else
             oldSX = StartX;
             oldSY = StartY;
             StartX = nextNeighborH.x;  %update nodes
             StartY = nextNeighborH.y;
             pointS(1,sIt) = convertedIndexS;
             sIt = sIt + 1;
        end
            
         if  intensityE > threshold
            mapPoints(convertedIndexE) = intensityE;
         else
           oldEX = EndX;
           oldEY = EndY;
           EndX = nextNeighbotT.x;
           EndY = nextNeighbotT.y;
           pointE(1,eIt) = convertedIndexE;
           eIt = eIt + 1;
         end
          sfoundPoints(1,it+1) =  convertedIndexS;
          efoundPoints(1,it+1) =  convertedIndexE;
     
           headOrTail = false;
           it = it +1;
           if abs(StartX- EndX)<= maskDim && ...
              abs(StartY-EndY)<= maskDim
               converged = 1;
             break  %converged
           end 

   
    end
   
   
    result.converged = converged;
    mergedPoints  = [pointS pointE];
    mergedPoints  = mergedPoints(mergedPoints > 0);
    mergedPointsSet = unique(mergedPoints,'rows','stable');
    path = zeros(1,size(mergedPointsSet,2));
    visitedPoints = zeros(1,size(mergedPointsSet,2));
    path(1,1) = pointS(1,1);
   % path(1,2) = pointE(1,1);
    %visitedPoints(1,1) = pointS(1,1);
    pathResult = finalPath(im,pointS(1,1),pointE(1,1),path,visitedPoints,2,1,mergedPointsSet,maskDim);
    pathResultPoints = pathResult.path;
    pathResultPoints  = flip(pathResultPoints(pathResultPoints > 0));
    pathResultPoints(1,size(pathResultPoints,2)) = pointE(1,1);
    pathResultPoints = [pointS(1,1) , pathResultPoints];  
    [xParts,yParts] = ind2sub(size(im), pathResultPoints);
    pointsMat = zeros(size(xParts,2),2);
    pointsMat(:,1)=xParts;
    pointsMat(:,2)=yParts;
    pointsMat =  unique(pointsMat,'rows','stable');
    result.points = pointsMat;
    %{    
     copyPathResults = pathResultPoints(1,:);
     pointsMat = zeros(size(xParts,2),2);
    
     pointsMat(:,1)=xParts;
     pointsMat(:,2)=yParts;
 
   % sortedPoints = unique(sortrows(mat),'rows');
%    result.points = mat;
    sampledPoints = zeros(14,2);
    currentPointX  = xParts(1);
    currentPointY  = yParts(1);
 %   currentIndex   =  pointS(1,1);
    sampledPoints(1,1) = currentPointX;
    sampledPoints(1,2) = currentPointY;
  %  skippedPoints  = zeros(1,size(pathResultPoints,2));
 %   skippedPoints(1,1) =  pointS(1,1);
    maxCords=sub2ind(size(im),size(im,1),size(im,2));
    copyPathResults(1,1) = maxCords;
    pointsMat(1,1) = size(im,1);
    pointsMat(1,2) = size(im,2);
   
    stepSize = 3;
    for c=1:14  % sample 14 points between start and end points from total detected points
        Offsets      = [-stepSize,0,stepSize];
         minInd = intmax;
         for xOffset = Offsets
           for yOffset = Offsets
               
               if ( xOffset == 0 ) && (yOffset == 0 )
                  continue
               end    
               newX    = (currentPointX+xOffset)-(floor(stepSize/2));
               newY    = (currentPointY+yOffset)-(floor(stepSize/2));
               [newX,newY] = adjustInvalidCordinates(newX,newY,stepSize);
               rangeResult = getInRangePoints(im,currentPointX,currentPointY,newX,newY,copyPathResults);
               if rangeResult.valid == false
                   continue
               end
                minInd =  rangeResult.minInd;
                I  = rangeResult.minIndLoc;
                copyPathResults(1,I)= maxCords;
                pointsMat(I,1) = size(im,1);
                pointsMat(I,2) = size(im,2);
                break
              % minVal = copyPathResults(1,I);
             %  if  any(skippedPoints(1,:) == minVal) 
              %       copyPathResults(1,I) = sub2ind(size(im),size(im,1),size(im,2));
              %       [minVal,I]  = min(abs(copyPathResults - ind));
              %        minVal = copyPathResults(1,I);
              % end
             %  if minVal < minInd
             %      minInd = minVal;
             %  end
           end
         end
         [nX,nY] = ind2sub(size(im),minInd);
         rangeResult = getInRangePoints(im,currentPointX,currentPointY,newX,newY,copyPathResults);
         if rangeResult.valid == true
             copyPathResults(1,rangeResult.Indices)= maxCords;  
             pointsMat(rangeResult.Indices,1) = size(im,1);
             pointsMat(rangeResult.Indices,2) = size(im,2);
         end

         %currentIndex  = minInd;
         [currentPointX,currentPointY] = ind2sub(size(im),minInd);
         sampledPoints(c+1,1) = currentPointX;
         sampledPoints(c+1,2) = currentPointY;
    end
    [EX,EY] = ind2sub(size(im),pointE(1,1));
    sampledPoints(16,1) = EX;
    sampledPoints(16,2) = EY;
    result.points = sampledPoints;
%}
 function [result]=getInRangePoints(im,currentPointX,currentPointY,nX,nY,pathResults)
     % sumMat     = ones(1,size(pathResults,2))*(size(im,1)+size(im,1));
      pointsMat  = zeros(size(pathResults,2),2);
      [xParts,yParts] =ind2sub(size(im),pathResults);
      pointsMat(:,1)=xParts;
      pointsMat(:,2)=yParts;
      filterPoints = pointsMat(pointsMat(:,1) >= min(currentPointX,nX),:);
      if  isempty(filterPoints) == 0
        filterPoints = filterPoints(filterPoints(:,1) <= max(currentPointX,nX),:);
      end
      if  isempty(filterPoints) == 0
         filterPoints = filterPoints(:,filterPoints(:,2) >= min(currentPointY,nY));
      end
      if  isempty(filterPoints) == 0 
         filterPoints = filterPoints(:,filterPoints(:,2) <= max(currentPointY,nY));
      end
      if isempty(filterPoints) == 1
          result.valid = false;
      else
          result.valid=true;
          [tf, loc] = ismember(filterPoints,pointsMat,'rows');
          result.Indices = loc;
          sumMat =  filterPoints(:,1) +  filterPoints(:,2);
          [minVal,I]  = min(abs(sumMat - (nX+nY)));
          chosenX   = filterPoints(I,1);
          chosenY   = filterPoints(I,2);
          [tf, loc] = ismember([chosenX chosenY],pointsMat,'rows');
          ind = sub2ind(size(im),chosenX,chosenY);
          result.minIndLoc = loc;
          result.minInd = ind;
      end
         
 
 function  [result]=finalPath(im,start,endP,path,visitedPoints,pIndex,vIndex,foundPoints,maskDim)
    
      
     Offsets      = [-maskDim,0,maskDim];
     neighbors    = zeros(1,8);
     [X,Y]        = ind2sub(size(im), start);
     nCount = 0;
     for xOffset = Offsets
        for yOffset = Offsets
             if ( xOffset == 0 ) && (yOffset == 0 )
                     continue
             end
             newX    = (X+xOffset)-(floor(maskDim/2));
             newY    = (Y+yOffset)-(floor(maskDim/2));
             [newX,newY] = adjustInvalidCordinates(newX,newY,maskDim);  
             neighborPoint = sub2ind(size(im), newX, newY);
                    
              if  any(visitedPoints(1,:) == neighborPoint)
                  continue
              end
            
             if neighborPoint == endP
                 result.valid = true;
                 result.visitedPoints = visitedPoints;
                 result.vIndex        = vIndex;
                 path(1,pIndex) =  start;  
                 pIndex = pIndex + 1;
                 result.path = path;
                 result.pIndex = pIndex;
                 return
             end
             
            
         
             if any(foundPoints(1,:) == neighborPoint)  
                 nCount = nCount +1;
                 neighbors(1,nCount) = neighborPoint; 
                 visitedPoints(1,vIndex) = neighborPoint;
                 vIndex = vIndex +1;  
             end
        end
     end
     
     if nCount== 0
         result.valid = false;
         visitedPoints(1,vIndex) = start;
         vIndex = vIndex +1;  
         result.pIndex = pIndex;
         result.path   = path;
         result.vIndex = vIndex;
         result.visitedPoints = visitedPoints;
         return
     end
     
    
     sumX=0;
     sumY=0;
     countValids = 0;
     result.valid = false;
     for i=1:nCount
         pathRes = finalPath(im,neighbors(1,i),endP,path,visitedPoints,pIndex,vIndex,foundPoints,maskDim);
         path =  pathRes.path;
         visitedPoints = pathRes.visitedPoints;
         pIndex = pathRes.pIndex;
         vIndex = pathRes.vIndex;        
         [nX,nY] = ind2sub(size(im), neighbors(1,i));  
         if pathRes.valid == true
             sumX = sumX + nX;
             sumY = sumY + nY;
             countValids = countValids +1;
             result.valid = true;
         end
     end

     visitedPoints(1,vIndex) = start;
     vIndex = vIndex +1;  
     result.vIndex = vIndex;
     result.visitedPoints = visitedPoints;
     if result.valid == true
         fX = round (sumX / countValids);
         fY = round (sumY / countValids);
         fPoint = sub2ind(size(im), fX, fY);
         path(1,pIndex) = fPoint;
         pIndex = pIndex + 1;
     end
     result.pIndex = pIndex;
     result.path   = path;


    
     
     
     %{  
    maxIteration =  size(mergedPointsSet,2);
    currentPoint  = pointS(1,1);
    visitedPoint = zeros(1,2*maxIteration); 
    resultPoints = zeros(1,2*maxIteration);
    neighbors    = zeros(1,maxIteration);
    Offsets = [-maskDim,0,maskDim];
    resultPoints(1,1) = currentPoint;
    visitedPoint(1,1) = currentPoint;
    it = 1;
    visitedIt = 1;
    targetPoint = pointE(1,1);  
    for c=1:2      
     while currentPoint ~= targetPoint  && it ~= maxIteration
      [X,Y] = ind2sub(size(im), currentPoint);
      i = 0;
      cX=0;
      cY=0;
      for xOffset = Offsets
        for yOffset = Offsets
         if ( xOffset == 0 ) && (yOffset == 0 )
                 continue
         end
         newX    = (X+xOffset)-(floor(maskDim/2));
         newY    = (Y+yOffset)-(floor(maskDim/2));
         [newX,newY] = adjustInvalidCordinates(newX,newY,maskDim);
         iPoint = sub2ind(size(im), newX, newY);
         if  any(visitedPoint(1,:) == iPoint)
             continue
         end
         if  any(mergedPointsSet(1,:) == iPoint)
            neighbors(1,i+1) = iPoint; 
            visitedIt = visitedIt + 1;
            visitedPoint(1,visitedIt) = iPoint;
          
            cX = cX + newX;
            cY = cY + newY;
            i = i+1;
         end
        end
      end
      if cX == 0 || cY == 0 
          break
      end
      averageX = round(cX / i);
      averageY = round(cY / i);
      convertedI = sub2ind(size(im), averageX, averageY);
      resultPoints(1,it+1) = convertedI;      
      currentPoint = convertedI;      
      it = it + 1 ;
      
     end  
      it = it + 1 ;
      resultPoints(1,it)  = pointE(1,1);  
      visitedIt = visitedIt + 1;
      visitedPoint(1,visitedIt) = pointE(1,1);
      currentPoint = pointE(1,1);
      targetPoint  = pointS(1,1); 
    end
    
    resultPoints  = resultPoints(resultPoints > 0);
    [xParts,yParts] = ind2sub(size(im), resultPoints);    
    mat = zeros(size(xParts,2),2);
    mat(:,1)=xParts;
    mat(:,2)=yParts;
    sortedPoints = unique(sortrows(mat),'rows');   
    result.points = sortedPoints; 
  %}
  
   

    
function[point] = findNeighbor(im,X,Y,maskDim,threshold,olderX,olderY,headOrTail,foundPoints,mapPoints)

mask = ones(maskDim, maskDim);
Offsets = [-maskDim,0,maskDim];

%if false
%    i = 0;
%    htOffsets = [-maskDim,0];
%    htNeighbors = ones(power(size(htOffsets),2));
%    for xOffset = htOffsets
%        for yOffset = htOffsets
%           i = i+1;
%           newX    = (X+xOffset);
%           newY    = (Y+yOffset); 
%           patch   = im(newX:newX+maskDim-1,newY:newY+maskDim-1);
%           intensity = sum(dot(patch,mask))/(maskDim*maskDim);  % convert it between 0 and 1   
%           htNeighbors(i) = intensity;
%        end
%    end
%    [val,Ind] = min(htNeighbors);
%    pointOne  = convertToPoint(X,Y,Ind,htOffsets);
%    point.x = pointOne.x+1;
%    point.y = pointOne.y+1;
%    point.intensity = val;
%    return
    
%end

neighbors = ones(power(size(Offsets),2));
weights = zeros(power(size(Offsets),2));
i = 0;
for xOffset = Offsets
   for yOffset = Offsets
     i =  i+1;
     if ( xOffset == 0 ) && (yOffset == 0 )
         continue
     end    
     newX    = (X+xOffset)-(floor(maskDim/2));
     newY    = (Y+yOffset)-(floor(maskDim/2));
     if newX == olderX && newY == olderY 
        continue
     end
     
     [newX,newY] = adjustInvalidCordinates(newX,newY,maskDim);
     patch   = im(newX:newX+maskDim-1,newY:newY+maskDim-1);
     weights(i) = nnz(patch <= threshold);
     intensity = sum(dot(patch,mask))/(maskDim*maskDim);  % convert it between 0 and 1   
     neighbors(i) = intensity;
   end
end

nCount=1;
while true
    [val,Ind] = min(neighbors);
     pointOne  = convertToPoint(X,Y,Ind,Offsets);
     iPoint = sub2ind(size(im), pointOne.x, pointOne.y);
     if  any(foundPoints(1,:) == iPoint) %||  isKey(mapPoints,iPoint)
        neighbors(Ind) = 1;
     else
         break
     end
     if nCount > size(neighbors)
         break
     end
     nCount = nCount + 1;
end
if false%~headOrTail 
    neighbors(Ind) = 1;
   [secondVal,secondInd] = min(neighbors);
   pointTwo  = convertToPoint(X,Y,secondInd,Offsets);
   distanceOne = sqrt(power(pointOne.x-olderX,2)+power(pointOne.y-olderY,2)); 
   distanceTwo = sqrt(power(pointTwo.x-olderX,2)+power(pointTwo.y-olderY,2));
   if distanceTwo > distanceOne 
       pointOne = pointTwo;
       val = secondVal;
   end

    
end


point.x = pointOne.x;
point.y = pointOne.y;
point.intensity =  im(pointOne.x,pointOne.y);




function[cpoint]=convertToPoint(X,Y,Ind,Offsets)
    chosenXOff = ceil(Ind / size(Offsets,2));        
    chosenYOff = mod(Ind,size(Offsets,2));
    if chosenYOff == 0 
        chosenYOff = size(Offsets,2);
    end;    
    cpoint.x = X+Offsets(chosenXOff);
    cpoint.y = Y+Offsets(chosenYOff);
    

function [point]=adjustInvalidPoint(iPoint,maskDim)
   point.x = iPoint.x;
   point.y = iPoint.y;
   if point.x <= 0 
          point.x = maskDim;
   end
   if point.y <= 0
          point.y = maskDim;
   end
   
function [x,y]=adjustInvalidCordinates(ix,iy,maskDim)
   x = ix;
   y = iy;
   if x <= 0 
          x = maskDim;
   end
   if y <= 0
         y = maskDim;
   end   
        
 






