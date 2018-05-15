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

%% build grain boundary model
stepSize = 1;
[gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID_0,x,y,stepSize);

%% draw grain boundary, and make handles



%%
stepSize = 2;

[boundaryTF,~,neighborID,tripleTF,~] = find_one_boundary_from_ID_matrix(ID_0);
tripleLookup = [x(tripleTF>0),y(tripleTF>0)];

ind = boundaryTF>0;
gbPoints = [ID_0(ind), neighborID(ind), x(ind), y(ind)];

% sort the two grain IDs corresponding to gb points
t = gbPoints(:,[1,2]);
t = sort(t,2);
gbPoints(:,[1,2]) = t;

uniquePair = unique(gbPoints(:,[1,2]),'rows');    % unique grain pairs

% reset variables
gb_dir = [];
gb_s_pt = [];
pt_s_gb = [];
pt_pos = [];

ptCount = 1;    % count number of control points
for igb = 1:length(uniquePair)
    inds = (gbPoints(:,1)==uniquePair(igb,1))&(gbPoints(:,2)==uniquePair(igb,2));
    segPts = gbPoints(inds,[3,4]);     % points of this grain boundary segment
    
    % [1] determine end points of grain boundary. determine and record direction
    if range(segPts(:,1))>=range(segPts(:,2))
        segPts = sortrows(segPts,1);
        gb_dir{igb} = 'horizontal';
    else
        segPts = sortrows(segPts,2);
        gb_dir{igb} = 'vertical';
    end
    
    % determine up to 5 keypoints of this grain boundary segment for fitting 
    inds = round(linspace(1,size(segPts,1), min(5, size(segPts,1))));
    keyPts = segPts(inds,:);
    
    for jj=1:size(keyPts,1)
        pt = keyPts(jj,:);  % current point considered
        
        % if first or last of keypoints, it could be a triple point
        if (jj==1)||(jj==size(keyPts,1))
            ind_triple = find(pdist2(pt,tripleLookup) <= sqrt(stepSize),1,'first');
            if ~isempty(ind_triple)
                pt = tripleLookup(ind_triple,:);
            end
        end
        
        % determine is the point was (e.g., a triple point that was) already used before
        try
            [~,loc] = ismember(pt,pt_pos,'rows');
        catch
            loc = 0;
        end

        if loc > 0
            ipt = loc;
            ptCount = ptCount - 1;  % to conteract the ++ at the end of loop
        else
            ipt = ptCount;
        end
        
        % [2] record this grain boundary's point id
        try
            gb_s_pt{igb} = [gb_s_pt{igb}, ipt];
        catch
            gb_s_pt{igb} = ipt;
        end
        
        % [3] record this 'pt'
        pt_pos(ipt,:) = pt;
        
        % [4] record this point's grain boundary id
        try
            pt_s_gb{ipt} = [pt_s_gb{ipt}, igb];
        catch
            pt_s_gb{ipt} = igb;
        end
        
        % [at the end of loop] increment ipt
        ptCount = ptCount + 1;

    end

end





