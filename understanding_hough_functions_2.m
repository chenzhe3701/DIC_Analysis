close all;
% generate a map
[x,y] = meshgrid(100:200,100:200);
map = zeros(size(x));

% make some lines
[pt, ind, indR, indC] = grids_covered_by_line(x,y,[150,100],[150,200]);
map(ind) = 1;

[pt, ind, indR, indC] = grids_covered_by_line(x,y,[120,100],[180,200]);
map(ind) = -1;

[pt, ind, indR, indC] = grids_covered_by_line(x,y,[100 160],[200 180]);
map(ind) = 1;

% show the lines
map = map .* rand(size(x));
figureHandle_1 = myplot(x,y,map);

% hough transform, houghpeaks the peaks, houghlines the line segments
[h,t,r] = hough(map);
[tg,rg] = meshgrid(t,r);
% figureHandle_2 = myplot(tg,rg,h);
figureHandle_2 = figure;
imshow(h,[],'XData',t,'YData',r,'colormap',parula); axis on;
xlabel('theta');ylabel('rho');

peaks = houghpeaks(h,32,'NHoodSize',[floor(size(h,1)/2)*2+1,3]);
set(0,'currentfigure',figureHandle_2); hold on;
for k = 1:size(peaks,1)
    xy = [tg(peaks(k,1),peaks(k,2)),rg(peaks(k,1),peaks(k,2))];
    plot3(xy(1),xy(2),max(h(:)),'s','LineWidth',1,'Color','k')
end



lines = houghlines(map,t,r,peaks);
% lines(k).point1/2 is in fact the index (index_c, index_r)
% show the extracted lines
set(0,'currentfigure',figureHandle_1); hold on;
for k = 1:length(lines)
    xy = [x(lines(k).point1(2),lines(k).point1(1)),y(lines(k).point1(2),lines(k).point1(1));...
        x(lines(k).point2(2),lines(k).point2(1)),y(lines(k).point2(2),lines(k).point2(1))];
    plot3(xy(:,1),xy(:,2),1*ones(size(xy,1),1),'LineWidth',2,'Color','green')
end