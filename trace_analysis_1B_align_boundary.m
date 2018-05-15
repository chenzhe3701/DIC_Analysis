% chenzhe, 2018-05-14
% In branch 'Ti7Al_B6'
% Organize this to look at grain boundary alignment

clear;
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

%%
th = quantile(exx(:),[0.005, 0.995]);
exx = mat_to_image(exx,th,'index');
%% build grain boundary model

indrs = 1:2000;
indcs = 1:2000;
exx_input = exx(indrs, indcs);

ID_input = ID(indrs, indcs);
x_input = X(indrs, indcs);
y_input = Y(indrs, indcs);
stepSize = y(2) - y(1);

[gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID_input,x_input,y_input,5);

%% draw grain boundary, and make handles
close all;
h = [];
H = [];
hline = [];
L = [];
G = [];
V = [];

figure;
imagesc([x_input(1),x_input(end)],[y_input(1),y_input(end)],exx_input);
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


%%
% stepSize = 2;
% 
% [boundaryTF,~,neighborID,tripleTF,~] = find_one_boundary_from_ID_matrix(ID_0);
% tripleLookup = [x(tripleTF>0),y(tripleTF>0)];
% 
% ind = boundaryTF>0;
% gbPoints = [ID_0(ind), neighborID(ind), x(ind), y(ind)];
% 
% % sort the two grain IDs corresponding to gb points
% t = gbPoints(:,[1,2]);
% t = sort(t,2);
% gbPoints(:,[1,2]) = t;
% 
% uniquePair = unique(gbPoints(:,[1,2]),'rows');    % unique grain pairs
% 
% % reset variables
% gb_dir = [];
% gb_s_pt = [];
% pt_s_gb = [];
% pt_pos = [];
% 
% ptCount = 1;    % count number of control points
% for igb = 1:length(uniquePair)
%     inds = (gbPoints(:,1)==uniquePair(igb,1))&(gbPoints(:,2)==uniquePair(igb,2));
%     segPts = gbPoints(inds,[3,4]);     % points of this grain boundary segment
%     
%     % [1] determine end points of grain boundary. determine and record direction
%     if range(segPts(:,1))>=range(segPts(:,2))
%         segPts = sortrows(segPts,1);
%         gb_dir{igb} = 'horizontal';
%     else
%         segPts = sortrows(segPts,2);
%         gb_dir{igb} = 'vertical';
%     end
%     
%     % determine up to 5 keypoints of this grain boundary segment for fitting 
%     inds = round(linspace(1,size(segPts,1), min(5, size(segPts,1))));
%     keyPts = segPts(inds,:);
%     
%     for jj=1:size(keyPts,1)
%         pt = keyPts(jj,:);  % current point considered
%         
%         % if first or last of keypoints, it could be a triple point
%         if (jj==1)||(jj==size(keyPts,1))
%             ind_triple = find(pdist2(pt,tripleLookup) <= sqrt(stepSize),1,'first');
%             if ~isempty(ind_triple)
%                 pt = tripleLookup(ind_triple,:);
%             end
%         end
%         
%         % determine is the point was (e.g., a triple point that was) already used before
%         try
%             [~,loc] = ismember(pt,pt_pos,'rows');
%         catch
%             loc = 0;
%         end
% 
%         if loc > 0
%             ipt = loc;
%             ptCount = ptCount - 1;  % to conteract the ++ at the end of loop
%         else
%             ipt = ptCount;
%         end
%         
%         % [2] record this grain boundary's point id
%         try
%             gb_s_pt{igb} = [gb_s_pt{igb}, ipt];
%         catch
%             gb_s_pt{igb} = ipt;
%         end
%         
%         % [3] record this 'pt'
%         pt_pos(ipt,:) = pt;
%         
%         % [4] record this point's grain boundary id
%         try
%             pt_s_gb{ipt} = [pt_s_gb{ipt}, igb];
%         catch
%             pt_s_gb{ipt} = igb;
%         end
%         
%         % [at the end of loop] increment ipt
%         ptCount = ptCount + 1;
% 
%     end
% 
% end





