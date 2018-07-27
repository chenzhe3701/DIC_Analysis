% script_test_onion_thin
%
% chenzhe, 2018-05-21


clear;
close all;
option = 4;

switch option
    case 0
        % option setting -1, load face image
        I = load('D:\p\m\DIC_Analysis\face.mat');
        I = I.I;
        I = double(I);
        
        J = ones(size(I));
    case 1
        % option setting - 2, draw two images
        I = zeros(255);
        J = zeros(255);
        
        figure; image(I);
        h = imellipse(gca, [25 25 220 188]);    % I is a ellipse
        I = h.createMask;
        I = double(I);
        close;
        
        % J is a polygon
        figure; image(J);
        h = impoly(gca, [76,24; 26,76; 26,186; 100,243; 161,217; 218,97; 200,80]);
        J = h.createMask;
        J = double(J);
        close;
        
    case 2
        % option setting - 2, draw two images
        I = zeros(255);
        J = zeros(255);
        
        figure; image(I);
        h = imellipse();
        % h = impoly();
        h.wait();
        
        I = h.createMask;
        I = double(I);
        close;
        
        % J is a polygon
        figure; image(J);
        % h = imellipse();
        h = impoly();
        h.wait();

        J = h.createMask;
        J = double(J);
        close;       
        
    case 3
        % run this to define 'I' using grain data
        ID = load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt\WE43_T6_C1_s_all_grain_1144_local_map.mat');
        ID = ID.data;
        ID = ID(5).ID_local;
        ID = ID==1144;
        A = double(ID);
        ind_start_input = find(A(:)==1,1,'first');
        A(:) = A(:).*[1:length(A(:))]';
        
        I = double(logical(ID));
%         I = A;
%         myplot(I);
        J = imwarp(I,projective2d([0.9 0.1 0.001;   0.07, 0.95, 0.001;   2, 5, 1]));
%         myplot(J);

    case 4
        % Run this to replace 'I' with grain data, then do some clean up,
        % because for a image with complex skeleton, the method does not
        % work well.  If we just simulate the grain with a simpler shape,
        % it works to some extend
        ID = load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt\WE43_T6_C1_s_all_grain_1144_local_map.mat');
        ID = ID.data;
        ID = ID(5).ID_local;
        ID = ID==1144;
        figure; image(zeros(size(ID)));
        h = impoly(gca, [257,14; 136,69; 14,298; 31,362; 253,503; 419,459; 428,203; 349,15]);
        I = h.createMask;
        I = double(I);
        close;
%         myplot(I);
        
        J = imwarp(I,projective2d([0.9 0.1 0.001;   0.07, 0.95, 0.001;   2, 5, 1]));
%         myplot(J);
        
        
end
%% plot the 2 images with skeleton labeled
[RGB1,R1,G1,B1] = color_canvas(size(I,1),size(I,2),1);
[RGB2,R2,G2,B2] = color_canvas(size(J,1),size(J,2),2);

I_skl = thin(I,inf);   % label the skeleton
II = I;
II(I_skl==1) = 0;
figure;image(cat(3,II,II,II).*RGB1);


J_skl = thin(J,inf);   % label the skeleton
JJ = J;
JJ(J_skl==1) = 0;
figure;image(cat(3,JJ,JJ,JJ).*RGB2);

%% perform thining
[pxl_M_1, ind_M_1, skl_M_1] = onion_thin(I);
[pxl_M_2, ind_M_2, skl_M_2] = onion_thin(J);

% plot the image using the sequence of layer traversal
figure; image(cat(3,R1(ind_M_1),G1(ind_M_1),B1(ind_M_1)));
figure; image(cat(3,R2(ind_M_2),G2(ind_M_2),B2(ind_M_2)));

% plot the skeleton point position during the traversal.
% Just to demonstrate that I have make the skeleton into single data point, unless all is skeleton
figure; imagesc(skl_M_1);
figure; imagesc(skl_M_2);

%% perform alignment
[xx,yy] = meshgrid(1:size(pxl_M_1,2), 1:size(pxl_M_1,1));

xq = linspace(1,xx(1,end),size(pxl_M_2,2));
yq = linspace(1,yy(end,1),size(pxl_M_2,1));

xq = repmat(xq, size(pxl_M_2,1), 1);
yq = repmat(yq', 1, size(pxl_M_2,2));

% if directly interp pixel value
% pxl_M_1_interp = interp2(xx,yy,pxl_M_1,xq,yq,'nearest');

% if interp index value
ind_M_1_interp = interp2(xx,yy,ind_M_1,xq,yq,'nearest');

I_interp = ones(size(J));
for iR = 1:size(ind_M_2,1)
    for iC = 1:size(ind_M_2,2)
        try
            I_interp(ind_M_2(iR,iC)) = ind_M_1_interp(iR,iC);
        catch
            I_interp(ind_M_2(iR,iC)) = 1;
        end
    end
end
figure;
% imagesc(img_interp);
image(cat(3,R1(I_interp),G1(I_interp),B1(I_interp)));






