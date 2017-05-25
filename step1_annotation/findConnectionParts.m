function [result]= findConnectionParts(im,StartX,StartY,EndX,EndY,maskDim)
    
    maxIteration = 100;
    pointS = zeros(maxIteration);
    pointE = zeros(maxIteration);
    it = 1;  % iterator
    while  it ~= maxIteration % to get sure not trapped  in infinity if they don't converge
        nextNeighborH = findNeighbor(im,StartX,StartY,maskDim); % starting from Head
        nextNeighbotT = findNeighbor(im,EndX,EndY,maskDim); % staring  from Tail
        StartX = nextNeighborH.x;  %update nodes
        StartY = nextNeighborH.y;
        EndX = nextNeighbotT.x;
        EndY = nextNeighbotT.y;
        pointS(it) = sub2ind(size(im), nextNeighborH.x, nextNeighborH.y);
        pointE(it) = sub2ind(size(im), nextNeighbotT.x, nextNeighbotT.y);
        it = it +1;
        if abs(nextNeighborH.x-nextNeighbotT.x)<= maskDim && ...
                abs(nextNeighborH.y-nextNeighbotT.y)<= maskDim
             break;  %converged

        end 
   
    end
    result.points  = [pointS;pointE];
    
function [point] = findNeighbor(im,X,Y,maskDim)

mask = ones(maskDim, maskDim);
Offsets = [-maskDim,0,maskDim];
neighbors = ones(size(Offsets)^2);
i = 1;
for xOffset = Offsets
   for yOffset = Offsets
     if ( xOffset == 0 ) && (yOffset == 0 )
         continue
     end    
     newX   = (X+xOffset)-(floor(maskDim/2));
     newY   = (Y+yOffset)-(floor(maskDim/2));
     patch  = im(newX:newX+maskDim,newY:newY+maskDim);  
     intensity = sum(dot(patch,mask))/(maskDim^2);  % convert it between 0 and 1   
     neighbors(i) = intensity;
     i =  i+1;
   end
end

[val,Ind] = min(neighbors);
chosenXOff = ceil(Ind / size(Offsets));         5 % 3 = 2  
chosenYOff = mod(Ind,size(Offsets));
if chosenYOff == 0 
    chosenYOff = size(Offsets);
end;    
point.x = X+chosenXOff;
point.y = Y+chosenYOff;
point.intensity = val;


 


