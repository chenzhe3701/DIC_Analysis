% chenzhe, 2018-03-20
%

% try to align (0) SEM-DIC during experiment (i.e., SEM with speckle), high
% reslution, but DIC data is only stepSize = 5, with correct [x,y] position
% (1) preTestSEM (stitched, medium resolution, missing one row) and
% (2) postTestSEM (stitched, medium resolution)

clear;
addChenFunction;
% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'boundaryTF');



dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

[filePreSEM,pathPreSEM] = uigetfile('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\preTestSEM\upper\img_full_res_.tif','select stitched SEM, preTest');
[filePreSEM_2,pathPreSEM_2] = uigetfile('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\preTestSEM\bottom row\img_full_res_.tif','select stitched SEM, preTest (lower part for WE43_T6_C1');
[filePostSEM,pathPostSEM] = uigetfile('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\postTestSEM\img_full_res_.tif','select stitched SEM, postTest');

%% load the images
Ipre = imread(fullfile(pathPreSEM,filePreSEM));
Ipre_2 = imread(fullfile(pathPreSEM_2,filePreSEM_2));
Ipost = imread(fullfile(pathPostSEM,filePostSEM));
load(fullfile(dicPath,dicFiles{end}),'exx','x','y','u','v');  

%% plot to find control points
figure;imshow(Ipre);
figure;imshow(Ipre_2);
figure;imshow(Ipost);
myplotm(exx,boundaryTF,'x',x,'y',y,'r',1);

%% control point position: [L1,R1; l1, r1; L2,R2,l2,r2]
cpPreSEM = [10003,2415; % L1
    51740,3860;     % R1
    13042,18400;    % l_1
    49000,17890;    % r_1
    ];
cpPreSEM_2 = [10720,564;      % L2
    46810, 1582;    % R2  
    10680,2500;     % l_2
    47030,3405;     % r_2
    ];

cpPostSEM = [5642,3360; % L1
    61540,6010;     % R1
    6905,36090;     % L2
    54600,37590;    % R2
    
    9821,27700;     % l_1
    57589,27160;    % r_1
    6727,39050;     % l_2
    54820,40370;    % r_2
];

% undeformed DIC
cpDIC = [2665, 2335;    % L1
    64370, 4415;        % R1
    3680, 33650;        % L2
    57050, 35070;       % R2
    
    7145,25760;     % l_1
    60360,24960;    % r_1
    3600,36450;     % l_2
    57320,37760;    % r_2
];

%% plot to check the selected control points
figure;imshow(Ipre);
hold on;
plot(cpPreSEM(:,1),cpPreSEM(:,2),'.r','markersize',64);

figure;imshow(Ipre_2);
hold on;
plot(cpPreSEM_2(:,1),cpPreSEM_2(:,2),'.r','markersize',64);

figure;imshow(Ipost);
hold on;
plot(cpPostSEM(:,1),cpPostSEM(:,2),'.r','markersize',64);

myplotm(exx,boundaryTF,'x',x,'y',y,'r',1);
hold on;
plot(cpDIC(:,1),cpDIC(:,2),'.r','markersize',64);

%% memory intensive, need to down-sample first

% (1) interp deformed_post_SEM to deformed_DIC. However, if plot using
% undeformed coordinates, this should make an undeformed SEM image using
% the deformed image.
cpDICp = zeros(size(cpDIC));
xp = x+u;
yp = y+v;
xp = inpaint_nans(xp);
yp = inpaint_nans(yp);    % there are lots of nans
for ii = 1:size(cpDIC,1)
   [ir,ic] = find((0==x-cpDIC(ii,1))&(0==y-cpDIC(ii,2)));
   cpDICp(ii,:) = [xp(ir,ic),yp(ir,ic)];
end
% can visualize to double check
% hold on;
% plot(cpDIC(:,1),cpDIC(:,2),'.r','markersize',64);
% plot(cpDICp(:,1),cpDICp(:,2),'.b','markersize',64);

tform = make_average_transform('projective',cpPostSEM,cpDICp);
% have to reduce size
ratio = round(size(Ipost,2)/size(x,2));
[xfrom,yfrom] = meshgrid(0:ratio:(size(Ipost,2)-1),0:ratio:(size(Ipost,1)-1));
Ipost = imresize(Ipost,[size(xfrom,1),size(xfrom,2)]);

Ipost_to_deformed_DIC = interp_data(xfrom,yfrom,Ipost,xp,yp,tform,'interp','nearest');
% IpostToUndeformed = Ipost_to_deformed_DIC; % if using pixel coordinates [x,y]
imwrite(Ipost_to_deformed_DIC,[saveDataPath,'IpostToUndeform.tif']);


% (2) upper part of pre_SEM to DIC/undeformed
tform = make_average_transform('affine',cpPreSEM,cpDIC([1,2,5,6],:));
% have to reduce size
ratio = round(size(Ipre,2)/size(x,2));
[xfrom,yfrom] = meshgrid(0:ratio:(size(Ipre,2)-1),0:ratio:(size(Ipre,1)-1));
Ipre = imresize(Ipre,[size(xfrom,1),size(xfrom,2)]);

Ipre_to_DIC = interp_data(xfrom,yfrom,Ipre,x,y,tform,'interp','nearest');
% IpostToUndeformed = Ipost_to_deformed_DIC; % if using pixel coordinates [x,y]
imwrite(Ipre_to_DIC,[saveDataPath,'IpreToDIC_1.tif']);


% (3) lower part of pre_SEM to DIC/undeformed
tform = make_average_transform('affine',cpPreSEM_2,cpDIC([3,4,7,8],:));
% have to reduce size
ratio = round(size(Ipre_2,2)/size(x,2));
[xfrom,yfrom] = meshgrid(0:ratio:(size(Ipre_2,2)-1),0:ratio:(size(Ipre_2,1)-1));
Ipre_2 = imresize(Ipre_2,[size(xfrom,1),size(xfrom,2)]);

Ipre_2_to_DIC = interp_data(xfrom,yfrom,Ipre_2,x,y,tform,'interp','nearest');
% IpostToUndeformed = Ipost_to_deformed_DIC; % if using pixel coordinates [x,y]
imwrite(Ipre_2_to_DIC,[saveDataPath,'IpreToDIC_2.tif']);


% (4) save exx as an image
clim = quantile(exx(:),[0.003,0.997]);
map = colormap;
close;
ind = interp1(clim, [1,size(map,1)], exx);
ind(isnan(ind))=1;
ind = round(ind);

R = reshape(map(ind(:),1),size(ind,1),size(ind,2));
G = reshape(map(ind(:),2),size(ind,1),size(ind,2));
B = reshape(map(ind(:),3),size(ind,1),size(ind,2));
RGB = cat(3,R,G,B);
imwrite(RGB,[saveDataPath,'img_of_exx.tif']);


