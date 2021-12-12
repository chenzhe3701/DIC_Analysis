% chenzhe 2016-4-12
% pt1 and pt2 defines a straight line
% x, y are matrcies defining the grid points
% pt are the (closest) points covered by this line
%
% chenzhe, 2017-09-01, fix the 'else' part bug.

function [pt,ind,indR,indC] = grids_covered_by_line(x,y,pt1,pt2)

% crop the local x and y matrix within this line region
xyLow = min([pt1;pt2],[],1);
xyHigh = max([pt1;pt2],[],1);
% index corresponding to the corner points
[~,ind1R] = min(abs(y(:,1)-xyLow(2)));
[~,ind1C] = min(abs(x(1,:)-xyLow(1)));
[~,ind2R] = min(abs(y(:,1)-xyHigh(2)));
[~,ind2C] = min(abs(x(1,:)-xyHigh(1)));
% local x,y data
xLocal = x(ind1R:ind2R,ind1C:ind2C);
yLocal = y(ind1R:ind2R,ind1C:ind2C);
% if line is not verticle
if pt1(1)~=pt2(1)
    k = (pt2(2)-pt1(2))/(pt2(1)-pt1(1));    % calculate k
    lineY = pt1(2) + k * (xLocal-pt1(1));   % based on a point's x-coordinate, it should have lineY as its y-coordinate if it's on the target line
    % now for each row or column, find the grid point closest to the target line 
    if abs(k)>=1
        [minValue,indC] = min(abs(yLocal-lineY),[],2); % for each row, find column index
        indR = (1:size(yLocal,1))';
    else
        [minValue,indR] = min(abs(yLocal-lineY),[],1); % for each column, find row index
        indC = (1:size(yLocal,2))';
    end
else
    lineY = yLocal; % if vertical, lineY is yLocal
    [minValue,indC] = min(abs(yLocal-lineY),[],2);
    indR = (1:size(yLocal,1))';
end
indR = reshape(indR,[],1);
indC = reshape(indC,[],1);
minValue = reshape(minValue,[],1);
% remove points outside of boarder, because it is the 'closest' to the line, but distance>half grid spacing, so it is not 'on' the line
[nR,nC] = size(yLocal);
uniqueValue = unique(y(:));
indRemove = (minValue>(uniqueValue(2)-uniqueValue(1))/2)&((indR==1)|(indR==nR)|(indC==1)|(indC==nC));
indR = indR(~indRemove);
indC = indC(~indRemove);
% calculate global
indR = indR + ind1R - 1;
indC = indC + ind1C - 1;

ind = sub2ind(size(x),indR,indC);
pt = [x(ind),y(ind)];





