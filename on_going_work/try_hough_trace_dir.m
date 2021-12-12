close all;

gData = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_GB_adjusted\WE43_T6_C1_s_all_grain_1144_local_map.mat';
data = load(gData);
data = data.data;
data = data(4);

X = data.x_local;
Y = data.y_local;
exx = data.exx_local;
exy = data.exy_local;
eyy = data.eyy_local;
cNumMap = data.clusterNumMapLocal;
ID = data.ID_local;
ID_current = data.ID_current;

myplot(cNumMap);
myplot(X,Y,exx);
%% perform analysis on binary cluster number map
close all;
map = cNumMap;
map(map~=1) = 0;
myplot(map);

[h,t,r] = hough(map);
[tgrid,rgrid] = meshgrid(t,r);
f = myplot(t,r,h);
axis normal;

%% try thining
close all;
thinMap = thin(map,inf);
myplot(thinMap);

[h,t,r] = hough(thinMap);
[tgrid,rgrid] = meshgrid(t,r);
myplot(t,r,h);
axis normal;

%% if perform on exx map
close all;
map = exx;
map(ID~=ID_current) = 0;
myplot(exx)

[h,t,r] = hough(map);
[tgrid,rgrid] = meshgrid(t,r);
myplot(t,r,h);
axis normal;



%% find peaks
nPeaks = 10;
peaks = houghpeaks(h,nPeaks);
hold on;
for k = 1:size(peaks,1)
    xy = [tgrid(peaks(k,1),peaks(k,2)),rgrid(peaks(k,1),peaks(k,2))];
    plot3(xy(1),xy(2),max(h(:)),'s','LineWidth',nPeaks+1-k,'Color','k')
end
