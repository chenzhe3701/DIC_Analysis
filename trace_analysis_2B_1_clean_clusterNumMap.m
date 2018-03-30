% for each grain, clean the cluster number map, eliminate small, un-connected region


clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
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

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
useParallel = 1;

%% select iEs to analyze cluster evolution
% [option 1] not use parallel
if ~useParallel
    for iE = iE_start:iE_stop
%         fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
        fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
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
            
            clusterNumMapLocal = one_pass_clean(clusterNumMapLocal);
            
            clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max) = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
            disp(['ID = ',num2str(ID_current)]);
        end
        
        save([saveDataPath,fName_c2t_result],'clusterNumMapCleaned','-append');
    end
end

%% [option 2] use parallel
if useParallel    
    for iE = iE_start:iE_stop
%         fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
        fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
        load([saveDataPath,fName_c2t_result],'clusterNumMap','stru');
        
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
                    
                    clusterNumMapLocal = one_pass_clean(clusterNumMapLocal);
                    
                    partMap{ii}(indR_min:indR_max, indC_min:indC_max) = partMap{ii}(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
                    disp(['ID = ',num2str(ID_current)]);
                end
            end
        end
        
        clusterNumMapCleaned = zeros(size(exx)); % new map of interest
        for ii=1:npool
            clusterNumMapCleaned = clusterNumMapCleaned + partMap{ii};
        end
        
        save([saveDataPath,fName_c2t_result],'clusterNumMapCleaned','-append');
    end    
end