close all;
I=imread('WE43_T6_C1_r0c0.tif');
I = imgaussfilt(I,3);
figure;imshow(I);


[~, threshold] = edge(I,'sobel');
fudgeFactor = .5;
BWs = edge(I,'sobel', threshold * fudgeFactor);
figure, imshow(BWs), title('binary gradient mask');

BWs = grow_boundary(BWs);
BWs = grow_boundary(BWs);
BWs = imfill(BWs,'holes');

figure, imshow(BWs), title('grown boundary');
%% 
BWs = one_pass_clean(BWs);
figure, imshow(BWs), title('cleaned boundary');

%%
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);
BWsdil = imdilate(BWs, [se90 se0]);
figure, imshow(BWsdil), title('dilated gradient mask');