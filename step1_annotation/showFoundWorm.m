function showFoundWorm(imOrg, wormFound, countWorm)
% show all the found worms in the image, annotated by segments on the parts
%
%
% see also rotflip2getMore, compact2full
% 
% Shu Kong
% 04/12/2015

%%
if ~exist('countWorm', 'var')
    countWorm = numel(wormFound);
end
%imshow(imOrg, 'border', 'tight');
imshow(imOrg); % ,'Border','tight'

for i = 1:countWorm
    line( wormFound{i}(2,1:end), wormFound{i}(1,1:end), 'linewidth', 1, 'color', 'r');
end