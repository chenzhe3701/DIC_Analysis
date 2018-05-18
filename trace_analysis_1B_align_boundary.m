% chenzhe, 2018-05-14
% In branch 'Ti7Al_B6'
% Organize this to look at grain boundary alignment

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
% resReduceRatio = 3;         % to save space, reduce map resolution
% grow_boundary_TF = 0;       % whether to grow boundary to make it thicker

% file name prefixes
f1 = 'Ti7Al_B6_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
% twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

%% modify exx, so that it is easier to be plotted by imagesc
th = quantile(exx(:),[0.005, 0.995]);
exx = mat_to_image(exx,th,'index');


%% crop an area of data
indrs = 2001:3600;
indcs = 3201:4800;
% indrs = 1:size(exx,1);
% indcs = 1:size(exx,2);
exx_input = exx(indrs, indcs);

ID_input = ID(indrs, indcs);
x_input = X(indrs, indcs);
y_input = Y(indrs, indcs);
stepSize = y_input(2) - y_input(1);

% show initial condition
gb = find_one_boundary_from_ID_matrix(ID_input);
myplot(exx_input);
myplot(ID_input);
myplot(exx_input, grow_boundary(gb));

%% build grain boundary model
resolution = 4096/120;
[gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID_input,x_input,y_input,resolution);

%% draw grain boundary, and make handles
close all;
h = []; H = []; hline = [];
L = []; G = []; V = [];

figure;
imagesc([x_input(1),x_input(end)],[y_input(1),y_input(end)],exx_input); % here to choose the background
a = gca;
hold on;

% (1) plot all the control points --> hangle: h{i}
for ii = 1:size(pt_pos,1)
   h{ii} = impoint(a, pt_pos(ii,1), pt_pos(ii,2));
   setColor(h{ii},'r');
end
% (2) for each boundary, group all its impoint handles --> gb_s_pt_group{j} = H{j} = {h{j1}, h{j2}, ... }
for jj = 1:length(gb_s_pt)
    H{jj} = h(gb_s_pt{jj}); % or, looks like this is the same { h{ gb_s_pt{jj} } } 
end
% (3) for each boundary, plot the line and record the handle: hline{j}
for jj = 1:length(gb_s_pt)
    hline{jj} = plot_spline_line(H{jj}, gb_dir{jj}, stepSize);    
end

% (4) for each point, group its related grain boundaries pt_s_gb_group{i} = L{i} = {hline{i1}, hline{i2}, ...}
% group this point's grain boundaries' points, pt_s_gb_s_pt_group{i} = G{i} = {H{i1}, H{i2}, ...}
% gropu this point's grain boundaries' direction, pt_s_gb_sdir{i} = V{i} = {gb_dir{i1}, gb_dir{i2}, ...}  
% addNewPositionCallback with L,G,V
for ii = 1:size(pt_pos,1)
    L{ii} = hline(pt_s_gb{ii});
    G{ii} = H(pt_s_gb{ii});
    V{ii} = gb_dir(pt_s_gb{ii});
    addNewPositionCallback(h{ii}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,stepSize) , L{ii}, G{ii}, V{ii}) );
end

%% update [pt_pos] --> maybe [gb_dir] --> even more, maybe [gb_s_pt] and [pt_s_gb]
for ii = 1:size(pt_pos,1)
    try
        pt_pos(ii,:) = h{ii}.getPosition;
    end
end

%% plot the mask
[mask,~] = plot_spline_mask(gb_dir, gb_s_pt, pt_pos, x_input, y_input);


%% After making some data, save the data to study how to walk one grain boundary to the target grain boundary.
save('boundary_mask.mat', 'mask', 'ID_input');
%% Can also save the boundary model
save('boundary_model.mat', 'gb_dir', 'gb_s_pt', 'pt_pos', 'pt_s_gb', 'tripleLookup','x_input','y_input','ID_input','exx_input','stepSize');

%% Load data
load('D:\p\m\DIC_Analysis\try_align_gb_data.mat','mask','ID_input');

gb_target = mask;
% d_target = city_block(gb_target);
% [FX,FY] = gradient(d_target);

% grow gb_target to guarantee grains are disconnected.  May need grow multiple times.   
gb_target = grow_boundary(gb_target);

% just show how ID_temp compares with old ID (but use grain boundary to show) 
myplot(ID_input, gb_target);


% find ID with boundary map. Temporarily make it negative, so it can be recognized if not matched.  
ID_temp = -find_ID_map_from_boundary_map(gb_target);

% just show how ID_temp compares with old ID (but use grain boundary to show)  
gb = find_one_boundary_from_ID_matrix(ID_input);
myplot(ID_temp, gb);


% change the id# in ID_temp to that in ID_input
ID_aligned = hungarian_assign_ID_map(ID_temp, ID_input);
gb_aligned = find_one_boundary_from_ID_matrix(ID_aligned);

% show the new, aligned ID map
myplot(ID_aligned, grow_boundary(gb_aligned));











