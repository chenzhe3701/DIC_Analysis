
% chenzhe, 2018-03-18
% Plan: read preEBSD and postEBSD SEM images, try different filters and see
% if I can reconstruct the grain boundary map from SEM images.

%% read image
close all;
I=imread('D:\WE43_T6_C1_insitu_compression\preEBSD_SEM\WE43_T6_C1_r0c0.tif');
I = I((1:2048)+1000, (1:2048)+1000);
figure;imshow(I);

Ip = imread('D:\WE43_T6_C1_insitu_compression\WE43_T6_C1_post\cropped\WE43_T6_C1_post_r0c0.tif');
Ip = Ip((1:2048)+1000, (1:2048)+1000);
figure;imshow(Ip);
%% read strain map
strainFile = 'D:\WE43_T6_C1_insitu_compression\byFov\r0c0\WE43_T6_C1_s5_r0c0.mat';
load(strainFile,'exx','exy','eyy','x','y','gamma');

%% (0) Gaussian smooth only
% (0.1) it is not recommended to use gaussian filter with fspecial, but why ?
m_sigma=3;
m_sz = 2*ceil(2*m_sigma) + 1;   % this is the default in imgaussianfilt
h = fspecial('gaussian', m_sz, m_sigma);
J = imfilter(I,h);
figure;
imshow(J);
title('using fspecial');

% (0.2) smooth image with gaussian filter
J = imgaussfilt(I,3);
figure;
imshow(J);
title('using imgaussfilt');

% to sum-up, they are almost the same. But maybe imgaussfilt uses the defautl sz=2*ceil(2*sigma)+1.
% Therefore, maybe for smoothing purpose, using imgaussfilt is more convenient.  

%% (1.1) edge function, with 'log' method
[gb,threshOut] = edge(I, 'log', [], 5);
figure;
imshow(gb);
title('edge() with log method');


%% (2) sobel method
J = imgaussfilt(I,3);
figure;
imshow(J);
title('gaussian smoothed using imgaussfilt');

[~, threshold] = edge(J,'sobel');
fudgeFactor = 1;
gb = edge(J,'sobel', threshold * fudgeFactor);
figure, imshow(gb), title('binary gradient mask');

%% post processing
gb = grow_boundary(gb);
gb = grow_boundary(gb);
figure, imshow(gb), title('grow boundary of the mask');


gb = imfill(gb,'holes');
figure, imshow(gb), title('fill holes with imfill');
%% 
% gb = one_pass_clean(gb);
% figure, imshow(gb), title('cleaned boundary');

%%
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
gb_dilated = imdilate(gb, [se90 se0]);
figure, imshow(gb_dilated), title('dilated gradient mask with imdilate');

se = strel('square',10);
gb_close = imclose(gb,se);
figure, imshow(gb_close), title('closed mask with imclose');





