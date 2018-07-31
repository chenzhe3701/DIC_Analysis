% chenzhe, 2018-05-14
% In branch 'Ti7Al_B6'
% Organize this to look at grain boundary alignment
%
% chenzhe, 2018-06-12. Previously, it is for development purpose.
% Next, we need to divide the map, and do the real work of clean-up.  

clear; clc;
addChenFunction;
dicPath = uigetdir('E:\Ti7Al_B6_insitu_tension\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('E:\Ti7Al_B6_insitu_tension\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'gID','gNeighbors','gNNeighbors','ID','ID_0','iThick','x','y','X','Y','exx','gExx');

gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 14;

% file name prefixes
f1 = 'Ti7Al_B6_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.

%%
% pre-deformation SEM image with grain boundary
[fileBoundaryImg, pathBoundaryImg] = uigetfile('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\ref imgs for gb\IpreToDIC_1.tif','select SEM img with grain boundary');
ImgGB = imread([pathBoundaryImg,fileBoundaryImg]);

%% crop an area of data
% indrs = 2001:3600;
% indcs = 3201:4800;
indrs = 1:size(exx,1);
indcs = 1:size(exx,2);
exx_input = exx(indrs, indcs);

ID_input = ID(indrs, indcs);
x_input = X(indrs, indcs);
y_input = Y(indrs, indcs);
stepSize = y_input(2) - y_input(1);

% show initial condition
gb = find_one_boundary_from_ID_matrix(ID_input);
myplot(exx_input);
myplot(ID_input);
myplot(exx_input, grow_boundary(grow_boundary(gb)));

%% choose strain level
loaded = load(fullfile(dicPath,'_11.mat'));
exx_input = loaded.exx;
th = quantile(exx_input(:),[0.005, 0.995]);
exx_input = mat_to_image(exx_input,th,'index');
exx_input = exx_input(indrs,indcs);

exy_input = loaded.exy;
th = quantile(exy_input(:),[0.005, 0.995]);
exy_input = mat_to_image(exy_input,th,'index');
exy_input = exy_input(indrs,indcs);

eyy_input = loaded.eyy;
th = quantile(eyy_input(:),[0.005, 0.995]);
eyy_input = mat_to_image(eyy_input,th,'index');
eyy_input = eyy_input(indrs,indcs);

% u_input = loaded.u; 
% u_input=u_input(indrs,indcs);
% 
% v_input = loaded.v; 
% v_input=v_input(indrs,indcs);
% 
% sigma_input = loaded.sigma; 
% sigma_input=sigma_input(indrs,indcs);
% 
% clear loaded;
% myplot(exy_input, (gb));
% myplot(eyy_input, (gb));
% myplot(u_input, (gb));
% myplot(v_input, (gb));
% myplot(sigma_input, (gb));

clear exx exy eyy;
%% build grain boundary model
% chenzhe, 2018-07-11, 'resolution' should be num_of_pixels/EBSD_step_size 
% Reason: gb smaller than this length is meaningless, because it was interp/extraped by an EBSD data point

img_resolution = 4096/360; % pixels / um
ebsd_resolution = 1/1; % # of EBSD_data_point / um
resolution = img_resolution / ebsd_resolution;

[gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID_input,x_input,y_input,resolution);
save([sampleName,'_boundary_model_temp.mat'], 'gb_dir', 'gb_s_pt', 'pt_pos', 'pt_s_gb', 'tripleLookup','stepSize');

%% Divide whole map into several AOIs
nPts = size(pt_pos,1);
nAOI = nPts/100;  % approximately 100 points per AOI
% determine [nr,nc] for AOI
rAOI = ceil(sqrt(nAOI/((x_input(end)-x_input(1))/(y_input(end)-y_input(1)))));
cAOI = ceil(rAOI * (x_input(end)-x_input(1))/(y_input(end)-y_input(1)));
% leave some extra space for AOI overlap
extraSize = (y(2)-y(1))*1000;
% find AOI window size
xWindowSize = (x_input(end)-x_input(1))/cAOI;
yWindowSize = (y_input(end)-y_input(1))/rAOI;
% define AOI{icell} = [x1,x2; y1,y2];
AOI = [];
for ir = 1:rAOI
   for ic = 1:cAOI
       AOI{ir,ic} = [floor(xWindowSize*(ic-1) - extraSize + x_input(1)), ceil(xWindowSize*ic + extraSize + x_input(1));...
           floor(yWindowSize*(ir-1) - extraSize + y_input(1)), ceil(yWindowSize*ir + extraSize + y_input(1))];
   end
end

%% Select an AOI by [ir,ic], use with modify_gb_model to update the grain boundary model
% Draw grain boundary and make handles, but use 'global index' for the points and grain boundaries
ir = 1;
ic = 1;

% find the x_input, y_input, exx_input (or exy_input, eyy_input) in this AOI  
[~,indrs_low] = min(abs(AOI{ir,ic}(2,1)-y_input(:,1)));
[~,indrs_high] = min(abs(AOI{ir,ic}(2,2)-y_input(:,1)));
[~,indcs_low] = min(abs(AOI{ir,ic}(1,1)-x_input(1,:)));
[~,indcs_high] = min(abs(AOI{ir,ic}(1,2)-x_input(1,:)));
indrs_AOI = indrs_low:indrs_high;
indcs_AOI = indcs_low:indcs_high;


%%%% select what to use as a background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

map_AOI = exx_input(indrs_AOI, indcs_AOI);

% map_AOI = mod(ID(indrs_AOI, indcs_AOI),8);

% map_AOI = ImgGB(indrs_AOI, indcs_AOI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


ID_AOI = ID(indrs_AOI, indcs_AOI);
x_AOI = X(indrs_AOI, indcs_AOI);
y_AOI = Y(indrs_AOI, indcs_AOI);
stepSize = y_AOI(2) - y_AOI(1);

% find all points in region, all points' boundaries
all_pts_ind = [];
all_gb_ind = [];
for ind_pts = 1:size(pt_pos,1)
   if (pt_pos(ind_pts,1)>AOI{ir,ic}(1,1))&&(pt_pos(ind_pts,1)<AOI{ir,ic}(1,2))&&(pt_pos(ind_pts,2)>AOI{ir,ic}(2,1))&&(pt_pos(ind_pts,2)<AOI{ir,ic}(2,2))
      all_pts_ind = [all_pts_ind, ind_pts];
      all_gb_ind = [all_gb_ind, pt_s_gb{ind_pts}]; 
   end
end
all_gb_ind = unique(all_gb_ind);
% find all additional points on the boundaries but maybe outside of the AOI   
for igb = all_gb_ind
   all_pts_ind = [all_pts_ind, gb_s_pt{igb}]; 
end
all_pts_ind = unique(all_pts_ind);

% draw and modify
close all;
h = []; H = []; hline = [];
L = []; G = []; V = [];

figure;
imagesc([x_AOI(1),x_AOI(end)],[y_AOI(1),y_AOI(end)],map_AOI); % here to choose the background
a = gca;
hold on;

% (1) plot all the control points --> hangle: h{i}
for ii = all_pts_ind %1:size(pt_pos,1)
   h{ii} = impoint(a, pt_pos(ii,1), pt_pos(ii,2));
   setColor(h{ii},'r');
end
% (2) for each boundary, group all its impoint handles --> gb_s_pt_group{j} = H{j} = {h{j1}, h{j2}, ... }
for jj = all_gb_ind %1:length(gb_s_pt)
    H{jj} = h(gb_s_pt{jj}); % or, looks like this is the same { h{ gb_s_pt{jj} } } 
    % ---------> also, maybe record the corresponding [point index].
end
% (3) for each boundary, plot the line and record the handle: hline{j}
for jj = all_gb_ind %1:length(gb_s_pt)
    hline{jj} = plot_spline_line(H{jj}, gb_dir{jj}, stepSize);  
end

% (4) for each point, group its related grain boundaries pt_s_gb_group{i} = L{i} = {hline{i1}, hline{i2}, ...}
% group this point's grain boundaries' points, pt_s_gb_s_pt_group{i} = G{i} = {H{i1}, H{i2}, ...}
% gropu this point's grain boundaries' direction, pt_s_gb_sdir{i} = V{i} = {gb_dir{i1}, gb_dir{i2}, ...}  
% addNewPositionCallback with L,G,V
for ii = all_pts_ind %1:size(pt_pos,1)
    pt_s_gb_in_aoi = intersect(pt_s_gb{ii}, all_gb_ind);
    L{ii} = hline(pt_s_gb_in_aoi);
    G{ii} = H(pt_s_gb_in_aoi);
    V{ii} = gb_dir(pt_s_gb_in_aoi);
    S{ii} = addNewPositionCallback(h{ii}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,p,ii,stepSize) , L{ii}, G{ii}, V{ii}) );
end
axis equal;

% can use this to find out close points
% for ii = all_pts_ind
%     ind = find((pdist2(pt_pos(ii,:),pt_pos(all_pts_ind,:))<80) & (pdist2(pt_pos(ii,:),pt_pos(all_pts_ind,:))>0));
%     if ~isempty(ind)
%         text(pt_pos(ii,1)+10,pt_pos(ii,2),'XXX')
%     end
% end

% can use this to save plot of an area
% set(gca,'xlim',[4500, 11000],'ylim',[6500,13000]);
% print('c:\users\zhechen\desktop\gray_post.tiff','-dtiff')

    
%% use modify_gb_model to modify and update [pt_pos, gb_dir, gb_s_pt, pt_s_gb] and all hangles and groups of hangles, etc   
% and can save temporarily
timeStr = datestr(now,'yyyymmdd_HHMM');
save([sampleName,'_boundary_model_temp_',timeStr,'.mat'], 'gb_dir', 'gb_s_pt', 'pt_pos', 'pt_s_gb', 'tripleLookup','stepSize');


%% plot the mask to check
[mask,~] = plot_spline_mask(gb_dir, gb_s_pt, pt_pos, x_input, y_input);
myplot(ID_input, grow_boundary(mask));  % check how it align with ID map. If not too different, that's fine 

% Save the plotted mask, and boundary model
save([sampleName,'_boundary_model.mat'], 'gb_dir', 'gb_s_pt', 'pt_pos', 'pt_s_gb', 'tripleLookup','x_input','y_input','ID_input','exx_input','stepSize','mask');


%%
gb_target = mask;

% grow gb_target to guarantee grains are disconnected.  May need grow multiple times.
gb_target = (grow_boundary(gb_target));

% just show how ID_temp compares with old ID (but use grain boundary to show) 
% myplot(ID_input, gb_target);
myplot(ID_input, grow_boundary(grow_boundary(gb_target)));

%% find ID with boundary map. Could make it [[[negative]]], so it can be recognized if not matched.  
ID_temp = find_ID_map_from_boundary_map(gb_target);

% just show how ID_temp compares with old ID (but use grain boundary to show)  
gb = find_one_boundary_from_ID_matrix(ID_input);
% myplot(ID_temp, (gb));
myplot(ID_temp, grow_boundary(grow_boundary(gb)));

%% change the id# in ID_temp to that in ID_input
[ID_aligned, ID_link_additional] = hungarian_assign_ID_map(ID_temp, ID_input, max(gID));

gb_aligned = find_one_boundary_from_ID_matrix(ID_aligned);

% show the new, aligned ID map
myplot(ID_aligned, grow_boundary(grow_boundary(grow_boundary(gb_aligned))));

% show the linked-ID of the newly-generated grains
maxInputID = max(ID_input(:));
ID_new_assign = ID_aligned;
ID_new_assign(ID_new_assign <= maxInputID) = 0;
for ii = 1:size(ID_link_additional,1)
    ID_new_assign(ID_new_assign==ID_link_additional(ii,1)) = ID_link_additional(ii,2);
end
myplot(ID_new_assign, grow_boundary(grow_boundary(gb_aligned)));
myplot(ID_input, grow_boundary(grow_boundary(gb)));


%% Appended the aligned_ID which is from the modeled grain boundary, to the gb model file
save([sampleName,'_boundary_model.mat'], 'ID_aligned', 'ID_link_additional', '-append');







