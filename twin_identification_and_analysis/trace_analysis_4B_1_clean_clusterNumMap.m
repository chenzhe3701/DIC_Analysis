% for each grain, clean the cluster number map, eliminate small, un-connected region
%
% chenzhe, 2018-09-15
% Based on 2B_1_clean_(), but use 'one_pass_fill' instead of 'one_pass_clean' to clean.   
% More modifications:  
% (1) Load data from 'on_the_fly', but save to 'cluster_result' at this step.   
% (2) Previously, in code 'B_2_threshold_cluster_()', we set clusters with size < 1000 with a negative cVol.
%     Maybe it's no longer necessary, so by default we don't run that clean up any more. 
%     However, we keep that option, and if we want to run, we just run it in this code.  
% (3) So we will run the code to determine cVol at this step.  
%
% chenzhe, 2018-12-17, based on 3B_1_clean_clusterNumMap
% try a newer method to clean, i.e., clean sub-clusters that are too close to the grain boundary  

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_2021_09','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','exx');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','exx');
end

% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));


% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 6;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
useParallel = 1;

threshold = 1000;
useThreshold = 0;

%% select iEs to analyze cluster evolution
% [option 1] not use parallel
if ~useParallel
    for iE = iE_start:iE_stop
        fName_source = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
        fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
        copyfile([saveDataPath,fName_source],[saveDataPath,fName_c2t_result],'f');

        load([saveDataPath,fName_c2t_result],'clusterNumMap','stru');
        
        clusterNumMapCleaned = zeros(size(exx)); % new map of interest
        for iS =1:length(stru)
            % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
            % iS = find(gIDwithTrace == 296); % for debugging.
            ID_current = gIDwithTrace(iS);
            
            
            ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
            indC_min = find(sum(ind_local, 1), 1, 'first');
            indC_max = find(sum(ind_local, 1), 1, 'last');
            indR_min = find(sum(ind_local, 2), 1, 'first');
            indR_max = find(sum(ind_local, 2), 1, 'last');
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
            
            clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
            clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
            
            clusterNumMapLocal = one_pass_fill_and_clean(clusterNumMapLocal, 0.001);
            
            clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max) = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
            disp(['ID = ',num2str(ID_current)]);
        end
        
        save(fullfile(saveDataPath,fName_c2t_result),'clusterNumMapCleaned','-append');
    end
end

%% [option 2] use parallel
if useParallel    
    for iE = iE_start:iE_stop
        fName_source = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
        fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
        copyfile(fullfile(saveDataPath,fName_source), fullfile(saveDataPath,fName_c2t_result),'f');
        
        load(fullfile(saveDataPath,fName_c2t_result),'clusterNumMap','stru');
        
        npool = 3;
        tN = 1:length(stru);
        a = npool - mod(length(tN),npool);
        tN = [zeros(1,a),tN];
        tN = reshape(tN,npool,[]);
        ipool = [];
        partMap = [];
        for ii = 1:npool
            ipool{ii} =  tN(ii,:);
            partMap{ii} = zeros(size(exx));
        end
        parfor ii = 1:npool
            for iS =ipool{ii}
                % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
                % iS = find(gIDwithTrace == 296); % for debugging.
                if iS > 0
                    ID_current = gIDwithTrace(iS);
                    
                    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
                    indC_min = find(sum(ind_local, 1), 1, 'first');
                    indC_max = find(sum(ind_local, 1), 1, 'last');
                    indR_min = find(sum(ind_local, 2), 1, 'first');
                    indR_max = find(sum(ind_local, 2), 1, 'last');
                    
                    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
                    
                    clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
                    clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
                    
                    clusterNumMapLocal = one_pass_fill_and_clean(clusterNumMapLocal, 0.001);
                    
                    partMap{ii}(indR_min:indR_max, indC_min:indC_max) = partMap{ii}(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
                    disp(['iPool: ',num2str(ii),', ID = ',num2str(ID_current)]);
                end
            end
        end
        
        clusterNumMapCleaned = zeros(size(exx)); % new map of interest
        for ii=1:npool
            clusterNumMapCleaned = clusterNumMapCleaned + partMap{ii};
        end
        
        save(fullfile(saveDataPath,fName_c2t_result),'clusterNumMapCleaned','-append');
    end    
end

%% This part summarizes the cluster volume.

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    
    load(fullfile(saveDataPath,fName_c2t_result),'clusterNumMap','stru','clusterNumMapCleaned');
    
    clusterNumMap = abs(clusterNumMap);
    clusterNumMapCleaned = abs(clusterNumMapCleaned);
    for iS = 1:length(stru)
        
        % initialize
        stru(iS).cVol = zeros(size(stru(iS).cLabel));   
        stru(iS).cVolCleaned = zeros(size(stru(iS).cLabel)); 
        
        for iCluster = 1:size(stru(iS).cLabel)
            cNum = stru(iS).cLabel(iCluster);
            
            ind = (ID==stru(iS).gID)&(abs(clusterNumMap)==cNum);
            sz = sum(ind(:));
            
            if useThreshold
                if sz > threshold
                    stru(iS).cVol(iCluster) = sz;
                else
                    stru(iS).cVol(iCluster) = -sz;
                    clusterNumMap(ind) = -cNum;
                    disp(['small cluster detected: ID=',num2str(stru(iS).gID),', cluster=',num2str(cNum),', size=',num2str(sz)]);
                end
            else
                stru(iS).cVol(iCluster) = sz;
            end
            
            ind = (ID==stru(iS).gID)&(abs(clusterNumMapCleaned)==cNum);
            sz = sum(ind(:));
            
            if useThreshold
                if sz > threshold
                    stru(iS).cVolCleaned(iCluster) = sz;
                else
                    stru(iS).cVolCleaned(iCluster) = -sz;
                    clusterNumMapCleaned(ind) = -cNum;
                    disp(['small cluster_cleaned detected: ID=',num2str(stru(iS).gID),', cluster=',num2str(cNum),', size=',num2str(sz)]);
                end
            else
                stru(iS).cVolCleaned(iCluster) = sz;
            end
        end
        if (rem(iS,100)==0)
            disp(['finished ',num2str(iS),' grains']);
        end
    end
    save(fullfile(saveDataPath,fName_c2t_result),'clusterNumMap','stru','clusterNumMapCleaned','-append');
end



