
% chenzhe, 2018-10-17
%
% find/label clusters that are close to gb, so they can be eliminated from
% analysis.

clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
% dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
% dicFiles = dir([dicPath,'\*.mat']);
% dicFiles = struct2cell(dicFiles);
% dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
end

% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

debugTF = 0;


%% Calculate the city block distance
cityDistMap = city_block(boundaryTF);

%% (0) Load cluster number data at all stops

cluster_number_maps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0
cluster_number_maps_cleaned = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    cluster_number_maps{iE} = clusterNumMap;
    cluster_number_maps_cleaned{iE} = clusterNumMapCleaned;
    
    % initialize/zero fields
    for iS =1:length(stru)
        stru(iS).cToGbDist = zeros(size(stru(iS).cLabel,1), 5);     % first element store gDia, then target_qts = [0.5, 0.8, 0.9, 0.95];
    end
    struCell{iE} = stru;
    cToGbDistMapCell{iE} = zeros(size(boundaryTF)); 
end

%% (1.1) match clusters and track volume evolution, fields related: cVol, cVolCleaned, preCluster, postCluster
hWaitbar = waitbar(0,'running each grain, each cluster, ...');

for iS =1:length(struCell{iE_start})
    % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
    
    ID_current = struCell{iE}(iS).gID;
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    cityDistLocal = cityDistMap(indR_min:indR_max, indC_min:indC_max);
    gDia = sqrt(sum(ID_local(:)==ID_current)/pi);
    
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
        
    for iE = iE_start:iE_stop
        cToGbDistMapLocal = zeros(size(ID_local));

        % cMap = cluster_number_maps{iE}(indR_min:indR_max, indC_min:indC_max);
        cMap_cleaned = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
        
        % cMap(ID_local~=ID_current) = 0;  % cluster number just this grain
        cMap_cleaned(ID_local~=ID_current) = 0;
        
        for iC=1:length(struCell{iE}(iS).cLabel)
            cNum = struCell{iE}(iS).cLabel(iC);
            ind = cMap_cleaned==cNum;
            ds = cityDistLocal(ind);
            qts = quantile(ds,[0.5, 0.8, 0.9, 0.95]);
            struCell{iE}(iS).cToGbDist(iC,:) = [gDia,qts];
            cToGbDistMapLocal(ind) = qts(end);
        end
        
        cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + cToGbDistMapLocal;
    end
    
    waitbar(iS/length(struCell{iE}), hWaitbar);
            
end

try
    close(hWaitbar);
end
warning('on','all');

%% (after 1.1) update to Save. Can disable to temporarily to prevent data loss.
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    cToGbDistMap = cToGbDistMapCell{iE};
    save([saveDataPath,fName_c2t_result],'stru','cToGbDistMap','-append');
end









