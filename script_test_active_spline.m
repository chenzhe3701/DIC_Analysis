clear;clc;
%%
I = imread('c:/Users/ZheChen/Desktop/untitled.tif');
I = I(101:500,201:600);
close all;
f = figure;

imagesc(I);
a = gca;
hold on;

h{1} = impoint(a,100,100);
h{2} = impoint(a,150,150);
h{3} = impoint(a,200,200);
h{4} = impoint(a,250,250);
h{5} = impoint(a,300,300);

stepSize = 5;

H{1} = {h{1},h{2},h{3}};    % each H{i} is for a grain boundary. Each h{i} is a control point.  
H{2} = {h{3},h{4},h{5}};


%% (test 2) test two grain boundaries, with a common point
for ii = 1:length(H)
    pos = H{ii}{1}.getPosition;
    xp = pos(1);
    for jj = 2:length(H{ii})
        pos = [pos; H{ii}{jj}.getPosition];
        xp = [xp, linspace(xp(end)+stepSize, pos(jj,1), 100)];
    end
    pp = csapi(pos(:,1),pos(:,2));
    yp = fnval(pp,xp);
    
    hline{ii} = plot(xp,yp,'--.');
end

for ii = 1:length(H)
    for jj = 1:length(H{ii})
        cell_of_line_handles = {hline{1},hline{2}};
        cell_of_cp_groups = {H{1}, H{2}};
        addNewPositionCallback(H{ii}{jj},@(p) cellfun(@(x,y) update_spline_line(p, x,y, stepSize) , cell_of_line_handles, cell_of_cp_groups ));
    end
end

%% (test 1) test a single grain boundary
pos = h{1}.getPosition;
xp = pos(1);
for ii = 2:length(h)
    pos = [pos; h{ii}.getPosition];
    xp = [xp, linspace(xp(end)+stepSize, pos(ii,1), 100)];
end
pp = csapi(pos(:,1),pos(:,2));


yp = fnval(pp,xp);
hline = plot(xp,yp,'--.k');

% add callback to update grain boundary line
for ii = 1:length(h)
    addNewPositionCallback(h{ii},@(p) update_spline_line(p, hline,h,stepSize));

end



