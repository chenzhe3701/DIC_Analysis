

% relabel twin, separated from previous trace_analysis_2().
% chenzhe, 2018-02-05

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

%% select iE to analyze
iE = 5;

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

% =========== match cluster with twin system again, if need to change parameter ============================
name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
load([saveDataPath,name_result_on_the_fly]);

scoreCF = 0.1;
if strcmpi(sampleName,'WE43_T6_C1')
    switch iE
        case 5
            scoreCF = 0.230;
        case 4
            scoreCF = 0.215;
        case 3
            scoreCF = 0.217;
        case 2
            scoreCF = 0.291;
    end
end
% scoreCF = 0.15;  % can manually modify
% criterion final, for modifying
% distCF = 0.035;
% sfCF = -1;
% shearTarget = 0.1289;
% shearCF = 0.05;
% costCF = 0.055;
msgbox('check scoreCF, iE, etc');

if ~exist('tNote','var')
    tNote.enable = [0 0];
    tNote.disable = [0 0];
end
%%
%     close all;
% Create a few maps to record the criterion.
twinMap = zeros(size(exx));
sfMap = zeros(size(exx));
disSimiMap = zeros(size(exx));
scoreMap = zeros(size(exx));

% twinMap_2 = zeros(size(exx));
% shearMap = zeros(size(exx));
% sfMap_2 = zeros(size(exx));
% costMap = zeros(size(exx));

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
    
    twinMapLocal = zeros(size(ID_local));          % local map to record twin_system_number
    sfMapLocal = zeros(size(ID_local));            % local map to record schmid_factor
    disSimiMapLocal = zeros(size(ID_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
    scoreMapLocal = zeros(size(ID_local));
    
%     twinMapLocal_2 = zeros(size(ID_local));          % local map to record twin_system_number
%     shearMapLocal = zeros(size(ID_local));            % local map to record schmid_factor
%     sfMapLocal_2 = zeros(size(ID_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
%     costMapLocal = zeros(size(ID_local));       % local map to record cost
    
    % ==== change cluster number into twin system number, or 0
    nCluster = length(stru(iS).cLabel);
    stru(iS).cEnable = zeros(nCluster,1);
    stru(iS).c2t = zeros(nCluster,1);
    for iCluster = 1:nCluster
        cNum = stru(iS).cLabel(iCluster);
        indClusterLocal = (clusterNumMapLocal==cNum);
        
        % ============== method-1 =============================
        pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
        %         pdistCS(pdistCS > distCF) = nan;                % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
        %         pdistCS(stru(iS).tSF < sfCF) = nan;             % [criterion-2] SF must > 0.35
        %         [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
        
        % score boundary y=kx+b passes [pdist2,sf] = [0, 0.15] and [0.05, 0.5]
        score = pdistCS * 7 + (0.5 - stru(iS).tSF);
        
        [m_score, ind_t] = nanmin(score,[],2);
        m_dist = pdistCS(ind_t);
        
        tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
        if ~isempty(tNote.enable) && ismember([ID_current,iCluster],tNote.enable,'rows')
            stru(iS).cEnable(iCluster) = 1;
        
        end
        if ~isempty(tNote.disable) && ismember([ID_current,iCluster],tNote.disable,'rows')
            stru(iS).cEnable(iCluster) = -1;
        end
 
        
        
%         % if this cluster in this grain is disabled, make it nan. If enabled, make score as 0.
%         if ~isempty(tNote.enable) && ismember([ID_current,iCluster],tNote.enable,'rows')
%             m_score = 0.00001;
%             m_dist = 0.00001;
%         end
%         if ~isempty(tNote.disable) && ismember([ID_current,iCluster],tNote.disable,'rows')
%             m_score = nan;
%             m_dist = nan;
%         end
        
%         if isnan(m_dist)
%             tsNum = [];
%         else
%             tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
%         end

        
        if tsNum > nss
            sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
            disSimiMapLocal(indClusterLocal) = m_dist;
            scoreMapLocal(indClusterLocal) = m_score;
            if (m_score < scoreCF) || (1 == stru(iS).cEnable(iCluster))
                twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map. For twinMap, assign if m_score < scoreCF
                scoreMapLocal(indClusterLocal) = m_score/10;
                stru(iS).c2t(iCluster) = tsNum;     % c2t is the identification label.  Cluster->Twin
            end
            if (-1 == stru(iS).cEnable(iCluster))
                twinMapLocal(indClusterLocal) = -tsNum;
                scoreMapLocal(indClusterLocal) = m_score*10;
                stru(iS).c2t(iCluster) = tsNum;
            end            
        end
        
%         % ================ method-2 ==========================
%         shearFit = stru(iS).cShear(iCluster,:);
%         shearFit(stru(iS).tSF < sfCF) = nan;            % [criterion-1] SF must > 0.35
%         shearFit(abs(shearFit-shearTarget) > shearCF) = nan;      % [criterion-2] shear_difference should be < 0.5
%         
%         costFit = stru(iS).cCost(iCluster,:);
%         costFit(costFit > costCF) = nan;                 % [criterion-3] costFit should be < 0.025.  Change to 0.035 empirically
%         
%         score = abs(shearFit-shearTarget) + 2*costFit;       % [criterion-4] the criterion is a combination of shearFit and sqrt(cost)
%         score(costFit*2+abs(shearFit-shearTarget) > scoreCF) = nan;      % [criterion-5] combination of shearDiff and sqrtCost should also be limited
%         
%         [m_score,ind_t] = min(score,[],2);
%         if isnan(m_score)
%             tsNum = [];
%         else
%             tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
%         end
%         if tsNum > nss
%             twinMapLocal_2(indClusterLocal) = tsNum;
%             sfMapLocal_2(indClusterLocal) = stru(iS).tSF(ind_t);
%             shearMapLocal(indClusterLocal)  = stru(iS).cShear(iCluster,ind_t);
%             costMapLocal(indClusterLocal) = stru(iS).cCost(iCluster,ind_t);
%         end
    end
    
    % copy identified twin system number to twinMap
    twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
    sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
    disSimiMap(indR_min:indR_max, indC_min:indC_max) = disSimiMap(indR_min:indR_max, indC_min:indC_max) + disSimiMapLocal;
    scoreMap(indR_min:indR_max, indC_min:indC_max) = scoreMap(indR_min:indR_max, indC_min:indC_max) + scoreMapLocal;
    
%     twinMap_2(indR_min:indR_max, indC_min:indC_max) = twinMap_2(indR_min:indR_max, indC_min:indC_max) + twinMapLocal_2;
%     shearMap(indR_min:indR_max, indC_min:indC_max) = shearMap(indR_min:indR_max, indC_min:indC_max) + shearMapLocal;
%     sfMap_2(indR_min:indR_max, indC_min:indC_max) = sfMap_2(indR_min:indR_max, indC_min:indC_max) + sfMapLocal_2;
%     costMap(indR_min:indR_max, indC_min:indC_max) = costMap(indR_min:indR_max, indC_min:indC_max) + costMapLocal;
    
    waitbar(iS/length(stru), hWaitbar);
    %     input('press to continue');
end

scoreMap(scoreMap==0)=nan;

try
    close(hWaitbar);
catch
end
dt = datetime;
save(['temp_result_s',num2str(iE),'_',num2str(dt.Year),'_',num2str(dt.Month),'_',num2str(dt.Day),'_',num2str(dt.Hour),'_',num2str(dt.Minute),'.mat'],'tNote','scoreCF');
save([saveDataPath,'twin_label_result_s',num2str(iE),'_',num2str(dt.Year),'_',num2str(dt.Month),'_',num2str(dt.Day),'_',num2str(dt.Hour),'_',num2str(dt.Minute),'.mat'],'iE','tNote','scoreCF');

%% adjust scale bar to select criterion.  run each of these individually as needed, and finally generate a twinMap. 
% (a)
[f,a,c,s,v]= myplotc(scoreMap,'x',X,'y',Y,'tf',boundaryTFB,'r',3);
%%
myplot(X, Y, twinMap,boundaryTFB);

%% save the result
name_result_modified = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_modified.mat'];
disp('start saving cluster_result_modified');
save([saveDataPath,name_result_modified],'clusterNumMap','twinMap','sfMap','disSimiMap','scoreMap','stru','tNote','scoreCF','sfCF');
% save([saveDataPath,name_result_modified],'stru','-append');
% save([saveDataPath,name_result_modified],'twinMap_2','shearMap','sfMap_2','costMap','-append');
disp('finished saving cluster_result_modified');


% temp code for plot and investigate the results
if 0
%     myplot(X,Y,exx,boundaryTFB);
%     myplot(X,Y,clusterNumMap,boundaryTFB);
%     
%     myplot(X,Y,disSimiMap,boundaryTFB);
%     myplot(X,Y,twinMap,boundaryTFB); caxis([18,24]);
%     myplot(X,Y,sfMap,boundaryTFB);
%     
%     myplot(X,Y,shearMap, boundaryTFB);
%     myplot(X,Y,abs(shearMap-0.1289), boundaryTFB); title('shear diff');
%     myplot(X,Y,twinMap_2, boundaryTFB); caxis([18,24]);
%     myplot(X,Y,sfMap_2, boundaryTFB);
%     myplot(X,Y,costMap, boundaryTFB);
end
