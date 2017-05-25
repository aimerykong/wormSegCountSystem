function mask = connPart4Body(mask, loc)
% Draw segments in the matrix/image to connect worm parts.
% Bresenham's line algorithm is used here for draw lines.
%
% Shu Kong
% 04/14/2015

%mask = mat2gray(patchTMP)*0;
%mask = zeros( size(mask,1), size(mask,2));

for PointPair = 1:size(loc,2)-1 % for all consecutive pairs of parts
    
    % get the two points in order, suppose x0<=x1
    x0 = double(loc(2,PointPair));
    x1 = double(loc(2,PointPair+1));
    y0 = double(loc(1,PointPair));
    y1 = double(loc(1,PointPair+1));
    if x0 > x1
        x0 = double(loc(2,PointPair+1));
        y0 = double(loc(1,PointPair+1));
        
        x1 = double(loc(2,PointPair));
        y1 = double(loc(1,PointPair));
    end
    
    
    %% Bresenham's line algorithm
    if x0 == x1 % special case -- vertical line
        for y = min([y0,y1]):max([y0,y1])
            mask(uint16(y), uint16(x0)) = 1;
        end        
    else
        deltax = x1 - x0;
        deltay = y1 - y0;
        error = 0;
        deltaerr = abs (deltay / deltax); % Assume deltax != 0 (line is not vertical),
        % note that this division needs to be done in a way that preserves the fractional part
        y = y0;
        for x = x0:x1
            mask(uint16(y), uint16(x)) = 1;
            error = error + deltaerr;
            while error >= 0.5
                mask(uint16(y), uint16(x)) = 1;
                y = y + sign(y1 - y0);
                error = error - 1.0;
            end
        end
    end    
end