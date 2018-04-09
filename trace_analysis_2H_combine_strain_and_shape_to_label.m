% relabel twin, separated from previous trace_analysis_2().
% chenzhe, 2018-02-05
% chenzhe, 2018-02-26, based on 2_label_twins_post().
%
% chenzhe, 2018-04-09 add note
% This code, based on manual clean results, labels twins and consider them
% as 'ground truth'.
% Fields added to 'stru'
%  stru(iS).dis = dissimilarity of cluster centroid to the twin strain of the best matching twin system 
%  stru(iS).sf = schmid factor of the best matching twin system
%  stru(iS).ts = twin system number of the best matching twin system 
%  stru(iS).trueTwin = confirmed (ground truth) twin system number 
%
% Also, generates a few maps and append to 'cluster_to_twin_result.mat'
% including 'strainScoreMap','shapeScoreMap','trueTwinMap'.


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

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

useCleanedMap = 0;      % in this code, cleaned or not is only related to the map drawn.  No training, or classification is used in this code.
% But, do pay attention to the image folder you selected.

%% select iE to analyze
for iE = iE_start:iE_stop
    %     iE = 3;
    
    % =========== match cluster with twin system again, if need to change parameter ============================
    % name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    % load([saveDataPath,name_result_on_the_fly]);
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    
    % can try to remove some fields
    stru = rmfield(stru,{'dis','sf','ts','trueTwin'});
    
    if useCleanedMap
        clusterNumMap = clusterNumMapCleaned;
    end
    % Create a few maps to record the criterion.
    cnnTwinMap = zeros(size(exx));
    %%%%%%%%% new map for TMS
    strainScoreMap = zeros(size(exx));
    shapeScoreMap = zeros(size(exx));
    trueTwinMap = zeros(size(exx));
    
    hWaitbar = waitbar(0,'Matching cluster with twin system ...');
    for iS =1:length(stru)
        
        ID_current = gIDwithTrace(iS);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        
        cnnTwinMapLocal = zeros(size(ID_local));
        %%%%%%%%%% new map for TMS
        strainScoreMapLocal = zeros(size(ID_local));
        shapeScoreMapLocal = zeros(size(ID_local));
        trueTwinMapLocal = zeros(size(ID_local));
        
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        
        %%%%%%%%%% for TMS summary new
        stru(iS).dis = zeros(nCluster,1);
        stru(iS).sf = zeros(nCluster,1);
        stru(iS).ts = zeros(nCluster,1);
        stru(iS).trueTwin = zeros(nCluster,1);
        
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            indClusterLocal = (clusterNumMapLocal==cNum);
            
            %%%%%%%%%% for TMS summary
            pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);
            [m_dist, ind_t] = nanmin(pdistCS,[],2);
            tsNum = stru(iS).tLabel(ind_t); 
            stru(iS).dis(iCluster) = m_dist;
            stru(iS).sf(iCluster) = stru(iS).tSF(ind_t);
            stru(iS).ts(iCluster) = tsNum;
            
            % c2t>0 means its either a twin, or satisfied the strainScore criterion but disabled.    
            % So, c2t>0 && cEnable>=0 means it is a real twin.
            if (stru(iS).c2t(iCluster)>0)&&(stru(iS).cEnable(iCluster)>=0)
                stru(iS).trueTwin(iCluster) = tsNum;
                trueTwinMapLocal(indClusterLocal) = tsNum;
            end
           	strainScoreMapLocal(indClusterLocal) = 7*stru(iS).dis(iCluster)-stru(iS).sf(iCluster);
            shapeScoreMapLocal(indClusterLocal) = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
          
        end
        
        % copy identified twin system number to twinMap
        strainScoreMap(indR_min:indR_max, indC_min:indC_max) = strainScoreMap(indR_min:indR_max, indC_min:indC_max) + strainScoreMapLocal;
        shapeScoreMap(indR_min:indR_max, indC_min:indC_max) = shapeScoreMap(indR_min:indR_max, indC_min:indC_max) + shapeScoreMapLocal;
        trueTwinMap(indR_min:indR_max, indC_min:indC_max) = trueTwinMap(indR_min:indR_max, indC_min:indC_max) + trueTwinMapLocal;
        
        waitbar(iS/length(stru), hWaitbar);
    end

    
    try
        close(hWaitbar);
    catch
    end
    
    disp('append cnn result to stru');
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    save([saveDataPath,fName_c2t_result],'stru','strainScoreMap','shapeScoreMap','trueTwinMap','-append')
end


