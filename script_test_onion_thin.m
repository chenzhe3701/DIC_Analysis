% script_test_onion_thin
%
% chenzhe, 2018-05-21



close all;
%% option setting -1, load face image
I = load('D:\p\m\DIC_Analysis\face.mat');
I = I.I;
I = double(I);

%% option setting - 2.1, draw first image
I = zeros(255,255);
figure;imagesc(I);
h = imellipse();
h = impoly();
I = h.createMask;
I = double(I).*reshape(1:length(I(:)),size(I,1),[]);
figure;imagesc(I); title('I');

J = ones(size(I));
%% option setting - 2.2, draw second image, instead of square
figure;imagesc(J);
h = imellipse();
J = h.createMask;
J = double(J).*reshape(1:length(J(:)),size(J,1),[]);
figure;imagesc(J); title('J');

%%
[pxl_M_1, ind_M_1, skl_M_1] = onion_thin(I);
[pxl_M_2, ind_M_2, skl_M_2] = onion_thin(J);

[x,y] = meshgrid(1:size(pxl_M_1,2), 1:size(pxl_M_1,1));

xq = linspace(1,x(1,end),size(pxl_M_2,2));
yq = linspace(1,y(end,1),size(pxl_M_2,1));

xq = repmat(xq, size(pxl_M_2,1), 1);
yq = repmat(yq', 1, size(pxl_M_2,2));

pxl_M_1_interp = interp2(x,y,pxl_M_1,xq,yq,'nearest');
ind_M_1_interp = interp2(x,y,ind_M_1,xq,yq,'nearest');

img_interp = zeros(size(J));
for iR = size(ind_M_2,1):-1:1
    for iC = size(ind_M_2,2):-1:1
        try
            img_interp(ind_M_2(iR,iC)) = pxl_M_1_interp(iR,iC);
        catch
            img_interp(ind_M_2(iR,iC)) = 1;
        end
    end
end
figure;imagesc(img_interp);