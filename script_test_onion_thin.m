% script_test_onion_thin
%
% chenzhe, 2018-05-21


clear; 
close all;
%% option setting -1, load face image
% I = load('D:\p\m\DIC_Analysis\face.mat');
% I = I.I;
% I = double(I);
% J = ones(size(I));

%% option setting - 2.1, draw first image
I = zeros(256,256);
figure;imagesc(I);
h = imellipse();
% h = impoly();
h.wait();
I = h.createMask;
I = double(I).*reshape(1:length(I(:)),size(I,1),[]);
I2 = thin(I,inf);
I3 = I; I3(I2)=0;
close;figure;imagesc(I3); title('I');

J = ones(size(I));
%% option setting - 2.2, draw second image, instead of square
figure;imagesc(J);
% h = imellipse();
h = impoly();
h.wait();
J = h.createMask;
J = double(J).*reshape(1:length(J(:)),size(J,1),[]);
J2 = thin(J,inf);
J3 = J; J3(J2)=0;
close;figure;imagesc(J); title('J');

%%
[pxl_M_1, ind_M_1, skl_M_1] = onion_thin(I);
[pxl_M_2, ind_M_2, skl_M_2] = onion_thin(J);


%%

figure; imagesc(pxl_M_1)
figure; imagesc(pxl_M_2)
figure; imagesc(skl_M_1)
figure; imagesc(skl_M_2)

%%
[x,y] = meshgrid(1:size(pxl_M_1,2), 1:size(pxl_M_1,1));

xq = linspace(1,x(1,end),size(pxl_M_2,2));
yq = linspace(1,y(end,1),size(pxl_M_2,1));

xq = repmat(xq, size(pxl_M_2,1), 1);
yq = repmat(yq', 1, size(pxl_M_2,2));

pxl_M_1_interp = interp2(x,y,pxl_M_1,xq,yq,'nearest');
ind_M_1_interp = interp2(x,y,ind_M_1,xq,yq,'nearest');

img_interp = zeros(size(J));
for iR = 1:size(ind_M_2,1)
    for iC = 1:size(ind_M_2,2)
        try
            img_interp(ind_M_2(iR,iC)) = pxl_M_1_interp(iR,iC);
        catch
            img_interp(ind_M_2(iR,iC)) = 1;
        end
    end
end
figure;imagesc(img_interp);


%% run this to replace 'I' with grain data
realGrain = 1;
if realGrain
    ID = load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt\WE43_T6_C1_s_all_grain_1144_local_map.mat');
    ID = ID.data;
    ID = ID(5).ID_local;
    ID = ID==1144;
    A = double(ID);
    ind_start_input = find(A(:)==1,1,'first');
end
A(:) = A(:).*[1:length(A(:))]';
I = A;
myplot(I);
J = imwarp(I,projective2d([0.9 0.1 0.001;   0.07, 0.95, 0.001;   2, 5, 1]));
myplot(J);





