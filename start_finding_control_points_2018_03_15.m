% chenzhe, 2018-03-14
% from post deformation img, using DIC data, construct a pre-deformation
% img.
% This might be one-time use.
%
% For WE43_T6_c1, The SEM map is unloaded, the SEM map at s7 is loaded.
% So, first, I should check if a projective transform can align the grain
% boundaries.
% If they align well, then I can use the post-deformation SEM + the DIC
% data to calculate a reference map at zero-strain, which can be used to
% try to generate a reference to clean EBSD data.

clear;
addChenFunction;
[f1,p1] = uigetfile('D:\WE43_T6_C1_insitu_compression\img_full_res_.tif','select the post deformation img');
[f2,p2] = uigetfile('D:\WE43_T6_C1_insitu_compression\stitched_DIC\_4.mat','select the DIC data at the last strain level');
img = imread(fullfile(p1,f1));
load(fullfile([p2,f2]),'exx','exy','eyy','x','y','u','v');

xData = 1:size(img,2);
yData = 1:size(img,1);

ratio = 5;
img = img(1:ratio:end,1:ratio:end);
xData = xData(1:ratio:end);
yData = yData(1:ratio:end);

%% use this to find out the control points

part = 'lr';
switch part
    case 'ul'
        ri = 1;
        rf = size(img,1)/4;
        ci = 1;
        cf = size(img,2)/4;
    case 'll'
        ri = size(img,1)/4*3;
        rf = size(img,1);
        ci = 1;
        cf = size(img,2)/4;
    case 'ur'
        ri = 1;
        rf = size(img,1)/4;
        ci = size(img,2)/4*3;
        cf = size(img,2);
    case 'lr'
        ri = size(img,1)/4*3;
        rf = size(img,1);
        ci = size(img,2)/4*3;
        cf = size(img,2);
end
xx = xData(ci:cf);
yy = yData(ri:rf);
imgp = img(ri:rf,ci:cf);

figure; 
imagesc(imgp,'xData',xx,'ydata',yy);
%%
myplot(x+u,y+v,exx);
%% control points
% compare to strain_7
cpPostSEM = [5625,3370;
    60730,4845;
    6940,36110;
    59150,32150;
];

cpDIC = [2625,2300;
    57470,4020;
    3745,34750;
    55810,31000;
];

% compare to strain_4.
% cpSEM = [6080, 4200;
% 60750, 5275;    
% 6905, 36080;
% 60940, 32600;
% ];
% 
% cpDIC = [3070, 3125;
%     62010, 3900;
% 3720, 33950;
% 62390, 30500;
% ];

t = maketform('projective',cpPostSEM,cpDIC);
% t = fitgeotrans(cpPostSEM,cpDIC,'projective');
% img2 = imwarp(img,t);

[xFrom,yFrom] = meshgrid(xData,yData);
xTo = x+u;  % interp [img] to position [xp yp], then if plot using [x y], it is directly the desired map!
yTo = y+v;
xTo = inpaint_nans(xTo);
yTo = inpaint_nans(yTo);
img2 = interp_data(xFrom,yFrom,img,xTo,yTo,t,'interp','nearest');

myplot(xTo,yTo,img2)
myplot(x,y,img2)
imwrite(img2,'WE43_T6_C1_reconstructed_undeformed_img.tif')






