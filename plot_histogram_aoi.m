% function [data_mean, data_std, mask] = hist_aoi(M,X(optional),Y(optional))
% X,Y = coordinate for the image
% M = data for the image
% 
% with a myplot() plotted map of strain fields/grain boundary overlays, ...
% manually crop an area, and plot histogram on an AOI on map
%
% chenzhe, 2018-01-08

function [data_mean, data_std, mask] = plot_histogram_aoi(M,varargin)

p = inputParser;
addRequired(p,'M');
addOptional(p,'X',repmat(1:size(M,2), size(M,1), 1));
addOptional(p,'Y',repmat([1:size(M,1)]', 1, size(M,2)));
parse(p,M,varargin{:})

M = p.Results.M;
X = p.Results.X;
Y = p.Results.Y;

limit_x_low = min(X(:));
limit_x_high = max(X(:));
limit_y_low = min(Y(:));
limit_y_high = max(Y(:));
f = figure;
hold on;
set(f,'position',[50,50,800,600]);

facealpha = 1;
imagesc([X(1),X(end)],[Y(1),Y(end)],M,'alphadata',facealpha);

clim = quantile(M(:),[0.005, 0.995]);
caxis(gca, [clim(1), clim(2)]);
c = colorbar;
axis equal;
set(gca,'xLim',[limit_x_low,limit_x_high],'yLim',[limit_y_low,limit_y_high]);

H = impoly;
H.wait;
mask = H.createMask;
delete(H);
data = M(mask);
figure;
h = histogram(data(:));

data_mean = nanmean(data(:));
data_std = nanstd(data(:));

legend(['mean:  ',num2str(data_mean),char(10),char(13),'std:     ',num2str(data_std)],'Location','best');
end

