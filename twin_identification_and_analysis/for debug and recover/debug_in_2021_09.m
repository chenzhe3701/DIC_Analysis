
clear;
clc;
addChenFunction;
% check cleaned cluster number map

% Load ID
% ID map
d = matfile('D:\WE43_T6_C1\Analysis_2021_09\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat');
ID = d.ID; 
boundary = find_one_boundary_from_ID_matrix(ID);
boundary = grow_boundary(boundary);

% ground truth data 'trueTwinMapCell'
d = matfile('D:\WE43_T6_C1\Analysis_2021_09\possibly useful data\WE43_T6_C1_trueTwinMapCell_published.mat');
trueTwinMapCell_ref = d.trueTwinMapCell;


% Initial identification result
d = matfile('D:\WE43_T6_C1\Analysis_2021_09\20211001_0027_twinMaps.mat');
initial_twin = d.twinMapCell;
struCell_initial = d.struCell;

% confirmed identification
d = matfile('D:\WE43_T6_C1\Analysis_2021_09\20211001_0212_relabeled_result.mat');
confirmed_twin = d.trueTwinMapCell;
struCell_confirmed = d.struCell;

% confirmed to twin variant, then cleaned, result
d = matfile('D:\WE43_T6_C1\Analysis_2021_09\WE43_T6_C1_new_variant_map.mat');
final_variant = d.variantMapCell;
final_cell = d.trueTwinMapCell;
struCell_final = d.struCell;

%% just show selected iE

saveDataPath = 'D:\WE43_T6_C1\Analysis_2021_09';

iE = 5;
% load cluster number map
fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
load(fullfile(saveDataPath,fName_c2t_result),'stru','clusterNumMap','clusterNumMapCleaned');

close all;
myplot(initial_twin{iE}, boundary); caxis([0 6]);
myplot(confirmed_twin{iE}, boundary); caxis([0 6]);
myplot(final_cell{iE}, boundary); caxis([18 24]);    %
myplot(trueTwinMapCell_ref{iE}, boundary); caxis([18 24]);


% check diff
diff = trueTwinMapCell_ref{iE} - final_cell{iE};
myplot(diff); caxis([-0.1 0.1])

sum(diff(:))
inds = diff~=0;
unique(ID(inds))


%%
clc;
inds = ID==1552;
myplot(ID==1552);

unique(initial_twin{iE}(inds))
unique(confirmed_twin{iE}(inds))

unique(final_cell{iE}(inds))
sum(final_cell{iE}(inds)==0)

unique(trueTwinMapCell_ref{iE}(inds))
sum(trueTwinMapCell_ref{iE}(inds)==0)


%% find local ind
ID_target = 193;
cluster_target = 4;
ind_local = ismember(ID, ID_target); %ismember(ID, [ID_current,ID_neighbor]);

[nR,nC] = size(ID);
indC_min = find(sum(ind_local, 1), 1, 'first');
indC_max = find(sum(ind_local, 1), 1, 'last');
indR_min = find(sum(ind_local, 2), 1, 'first');
indR_max = find(sum(ind_local, 2), 1, 'last');
% increase area size by 1 pixel. 2021-09-29
indC_min = max(indC_min-1, 1);
indC_max = min(indC_max+1, nC);
indR_min = max(indR_min-1, 1);
indR_max = min(indR_max+1, nR);

ID_local = ID(indR_min:indR_max, indC_min:indC_max);
cMapLocal = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max);
myplot(cMapLocal)

% this grain, and this cluster
ind = ismember(cMapLocal, cluster_target) & ismember(ID_local, ID_target);
% confirmed twin
temp = confirmed_twin{iE}(indR_min:indR_max, indC_min:indC_max);
myplot(temp);
temp(~ind) = 0;
myplot(temp);
% final twin
temp = final_variant{iE}(indR_min:indR_max, indC_min:indC_max);
myplot(temp);
temp(~ind) = 0;
myplot(temp);
% ground truth twin
temp = trueTwinMapCell_ref{iE}(indR_min:indR_max, indC_min:indC_max);
myplot(temp);
temp(~ind) = 0;
myplot(temp);

try
    iS = find(arrayfun(@(x) x.gID==ID_target, struCell_initial{2}) )
catch
    iS = find(arrayfun(@(x) x.gID==ID_target, struCell{2}) )
end















