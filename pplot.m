% function [] = pplot(x_in,y_in,z_in)
%
% input can be vectors. Reconstruct mesh to plot a map.
% chenzhe, 2018-01-08

function [] = pplot(x_in,y_in,z_in)
% find range
xs = sort(unique(x_in(:)));
xStep = xs(2)-xs(1);
ys = sort(unique(y_in(:)));
yStep = ys(2)-ys(1);
xmin = min(x_in(:));
xmax = max(x_in(:));
ymin = min(y_in(:));
ymax = max(y_in(:));
% make grid
[xMesh,yMesh] = meshgrid(xmin:xStep:xmax,ymin:yStep:ymax);
[nR,nC] = size(xMesh);
subC = (x_in-xmin)/xStep+1;
subR = (y_in-ymin)/yStep+1;
ind = sub2ind([nR,nC],subR,subC);

zValue = nan*zeros(nR,nC);
zValue(ind) = z_in;

myplot(xMesh,yMesh,zValue);