function wormSetMore = compact2full(wormSetMore)
% get the rotated and flipped image accordingly based on the compact
% storage.
%
%
%
% see also rotflip2getMore, showFoundWorm
%
% Shu Kong
% 04/13/2015

%% generate more annotated image by rotating the original image with its worm labels
count = 2;
imOrg = wormSetMore{1}.im;
for i = 2:4
    wormSetMore{count}.im = eval(wormSetMore{count}.im);
    count = count + 1;
end

imOrg = eval(wormSetMore{count}.im); 
wormSetMore{count}.im = imOrg;
count = count + 1;
for i = 2:4
    wormSetMore{count}.im = eval(wormSetMore{count}.im);
    count = count + 1;
end









