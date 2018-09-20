
% relabel twin, separated from previous trace_analysis_2().
% chenzhe, 2018-02-05
% chenzhe, 2018-02-26, based on 2_label_twins_post().

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% DIC data, sometimes useful
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx','gPhi1','gPhi','gPhi2','gCenterX','gCenterY');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx','gPhi1','gPhi','gPhi2','gCenterX','gCenterY');
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

saveTF = 1;
scoreForDisabled = 0;   % arbitrary, but depend on data distribution.  Assign to score for disabled clusters. 
scoreForEnabled = 1;

%% select iE to analyze
iE = 3;

% strain data
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile)
clear('exy_corrected');
load(strainFile,'exx','exy','eyy','sigma','exy_corrected');     % if 'exy_corrected' does not exist, this does not give error, rather, just warning. % ----------------------------------------------------------------------------------
if exist('exy_corrected','var')&&(1==exy_corrected)
    disp('================= exy already corrected ! ========================');
    exy_corrected = 1;
else
    disp('================= exy being corrected here ! =======================');
    exy = -exy;
    exy_corrected = 1;
end
% remove bad data points
exx(sigma==-1) = nan;
exy(sigma==-1) = nan;
eyy(sigma==-1) = nan;
qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
exx(ind_outlier) = nan;
exy(ind_outlier) = nan;
eyy(ind_outlier) = nan;
% end strain data    
    
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

% =========== match cluster with twin system again, if need to change parameter ============================
% name_source = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat']    % data from previous step.
name_source = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat']    % data from previous step.
load([saveDataPath,name_source]);

% Threshold values.
scoreCF = 0;   % to start with. After adjusting threshold, disable this line.
if strcmpi(sampleName,'WE43_T6_C1')
    switch iE
        case 5
            scoreCF = 0.964;    % -0.243;
        case 4
            scoreCF = 0.964;    % -0.218;
        case 3
            scoreCF = 0.964;    % -0.218;
        case 2
            scoreCF = 0.964;    % -0.218;
    end
end


% criterion final, for modifying
% distCF = 0.035;
% sfCF = -1;
% shearTarget = 0.1289;
% shearCF = 0.05;
% costCF = 0.055;
msgbox('check source file, scoreCF, iE, etc');

if ~exist('tNote','var')
    tNote.enable = [0 0];
    tNote.disable = [0 0];
end
%%
%     close all;
% Create a few maps to record the criterion.
twinMap = zeros(size(exx));
sfMap = zeros(size(exx));
mDistMap = zeros(size(exx));
scoreMap = zeros(size(exx));

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
    mDistMapLocal = zeros(size(ID_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
    scoreMapLocal = zeros(size(ID_local));
    
    % ==== change cluster number into twin system number, or 0
    nCluster = length(stru(iS).cLabel);
    stru(iS).cEnable = zeros(nCluster,1);
    stru(iS).c2t = zeros(nCluster,1);
    for iCluster = 1:nCluster
        % if cVolCleaned==0, move on to the next cluster. The following 3 lines is east to comment out if necessary. 
        if (isfield(stru,'cVolCleaned'))&&(stru(iS).cVolCleaned(iCluster)<=0)
            continue;
        end
        
        cNum = stru(iS).cLabel(iCluster);
        indClusterLocal = (clusterNumMapLocal==cNum);
        
        % ============== method-1 =============================
        pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
        %         pdistCS(pdistCS > distCF) = nan;                % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
        %         pdistCS(stru(iS).tSF < sfCF) = nan;             % [criterion-2] SF must > 0.35
        %         [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
        
        % A boundary y=kx+b passes [pdist2,sf] = [0, 0.15] and [0.05, 0.5], then map into [0, 1] to define the phi-classifier 
        % This is the straight-forward definition.
        score = pdistCS * 7 - stru(iS).tSF;
        
        [m_dist, ind_t] = nanmin(pdistCS,[],2);     % Select the best-match twin system. 
        m_score = score(ind_t);
        % The other way is:
        % [m_score, ind_t] = nanmin(score,[],2);
        % m_dist = pdistCS(ind_t);

        % Then convert m_socre to phi_classifier (phi_score)
        C = 7;  % C is slope
        H = 1;  % H is the truncate position for mDist 
        % With this, phi_classifier = (SF - 7*mDist)/8 + 15/16 
        phi_score = stru(iS).tSF(ind_t)/(1+C*H) -C/(1+C*H)*min(H,m_dist) + (1+2*C*H)/2/(1+C*H);
        
        tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
        if ~isempty(tNote.enable) && ismember([ID_current,iCluster],tNote.enable,'rows')
            stru(iS).cEnable(iCluster) = 1;
            
        end
        if ~isempty(tNote.disable) && ismember([ID_current,iCluster],tNote.disable,'rows') 
            if (phi_score > scoreCF)
                stru(iS).cEnable(iCluster) = -1;        % to qualify to 'disable', it has to be identified as twin if were not disabled (phi_score < scoreCF)
            end
        end
        
        
        if tsNum > nss
            sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
            mDistMapLocal(indClusterLocal) = m_dist;
            scoreMapLocal(indClusterLocal) = phi_score;
            if ((phi_score > scoreCF)&&(-1 ~= stru(iS).cEnable(iCluster))) || (1 == stru(iS).cEnable(iCluster))

                if phi_score > scoreCF
                    % (1) condition: if phi_score < scoreCF, then no need to use 'Enable', so set the .cEnable(iCluster) back to 0
                    stru(iS).cEnable(iCluster) = 0;
                    scoreMapLocal(indClusterLocal) = phi_score + 1;
                else
                    % (2) condition: else, it needs to be enabled
                    scoreMapLocal(indClusterLocal) = scoreForEnabled;
                end
                twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map. For twinMap, assign if phi_score < scoreCF
                stru(iS).c2t(iCluster) = tsNum;     % c2t is the identification label.  Cluster->Twin
            end
            if (-1 == stru(iS).cEnable(iCluster))
                % (3) condition: if disabled
                twinMapLocal(indClusterLocal) = -tsNum;
                scoreMapLocal(indClusterLocal) = scoreForDisabled;
                stru(iS).c2t(iCluster) = tsNum;
            end
        end
        
    end
    
    % copy identified twin system number to twinMap
    twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
    sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
    mDistMap(indR_min:indR_max, indC_min:indC_max) = mDistMap(indR_min:indR_max, indC_min:indC_max) + mDistMapLocal;
    scoreMap(indR_min:indR_max, indC_min:indC_max) = scoreMap(indR_min:indR_max, indC_min:indC_max) + scoreMapLocal;
    
    try
        waitbar(iS/length(stru), hWaitbar);
    catch
    end
end

scoreMap(scoreMap==0)=nan;

try
    close(hWaitbar);
catch
end

if saveTF
    timeStr = datestr(now,'yyyymmdd_HHMM');
    save(['temp_tNote_s',num2str(iE),'_',timeStr,'.mat'],'tNote','scoreCF');
    save([saveDataPath,'twin_label_tNote_s',num2str(iE),'_',timeStr,'.mat'],'iE','tNote','scoreCF');
end
%%
[~,ae,~] = myplotm(exx,'x',X,'y',Y,'tf',boundaryTFB,'r',2);
colormapA('parula');
% 
[~,at,~] = myplot(X, Y, twinMap,boundaryTFB); caxis([0 24]);
% %%
% myplot(X, Y, grow_boundary((shrink_boundary((twinMap>0)))),boundaryTFB);
% %%
% myplot(X, Y, grow_boundary(grow_boundary(shrink_boundary(shrink_boundary((twinMap>0))))),boundaryTFB);
%% adjust scale bar to select criterion.  run each of these individually as needed, and finally generate a twinMap.
% [f,a,c,s,v]= myplotc(scoreMap,'x',X,'y',Y,'tf',boundaryTFB,'r',2);
[f,a,c,s,v]= myplotc_high(scoreMap,'x',X,'y',Y,'tf',boundaryTFB,'r',2);     % This is used to select value higher than threshold
% colormapA;
cl = caxis;
caxis([scoreCF,max(scoreCF+0.1*(cl(2)-cl(1)),cl(2))])

%% try this for disable using twinMap
[f,a,c,s,v]= myplotc_low(twinMap,'x',X,'y',Y,'tf',boundaryTFB,'r',2);
caxis([-20,10]);

%% get ID from map, then plot unit cell to check trace, one at a time
ids = find_ID_on_map(X,Y,ID,f,a);
ind = find(gID==ids(1));
hcp_cell('euler',[gPhi1(ind),gPhi(ind),gPhi2(ind)], 'ss', 25:30, 'stress', [-1 0 0; 0 0 0; 0 0 0]);

%% label_grain of interest on map to observe
label_map_with_ID(X,Y,ID,gcf, 102);

%% Use tNote to modify
hWaitbar = waitbar(0,'Re-matching cluster with twin system ...');
for iS =1:length(stru)
    
    ID_current = gIDwithTrace(iS);
    
    % if need to modify this grain
    if (~isempty(tNote.enable) && ismember(ID_current,tNote.enable(:,1))) || (~isempty(tNote.disable) && ismember(ID_current,tNote.disable(:,1)))
        
        % reset on global map the region of this grain 
        ind_global = (ID==ID_current);
        twinMap(ind_global) = 0;
        sfMap(ind_global) = 0;
        mDistMap(ind_global) = 0;
        scoreMap(ind_global) = 0;
        
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        
        % to modify, first copy the old values.
        twinMapLocal = twinMap(indR_min:indR_max, indC_min:indC_max);
        sfMapLocal = sfMap(indR_min:indR_max, indC_min:indC_max);
        mDistMapLocal = mDistMap(indR_min:indR_max, indC_min:indC_max);
        scoreMapLocal = scoreMap(indR_min:indR_max, indC_min:indC_max);
        
        %     twinMapLocal_2 = zeros(size(ID_local));          % local map to record twin_system_number
        %     shearMapLocal = zeros(size(ID_local));            % local map to record schmid_factor
        %     sfMapLocal_2 = zeros(size(ID_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
        %     costMapLocal = zeros(size(ID_local));       % local map to record cost
        
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        stru(iS).cEnable = zeros(nCluster,1);
        stru(iS).c2t = zeros(nCluster,1);
        for iCluster = 1:nCluster
            % if cVolCleaned==0, move on to the next cluster. The following 3 lines is east to comment out if necessary.
            if (isfield(stru,'cVolCleaned'))&&(stru(iS).cVolCleaned(iCluster)<=0)
                continue;
            end
        
            cNum = stru(iS).cLabel(iCluster);
            indClusterLocal = (clusterNumMapLocal==cNum);
            
            % ============== method-1 =============================
            pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
            %         pdistCS(pdistCS > distCF) = nan;                % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
            %         pdistCS(stru(iS).tSF < sfCF) = nan;             % [criterion-2] SF must > 0.35
            %         [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
            
            % A boundary y=kx+b passes [pdist2,sf] = [0, 0.15] and [0.05, 0.5], then map into [0, 1] to define the phi-classifier
            % This is the straight-forward definition.
            score = pdistCS * 7 - stru(iS).tSF;
            
            [m_dist, ind_t] = nanmin(pdistCS,[],2);     % Select the best-match twin system.
            m_score = score(ind_t);
            % The other way is:
            % [m_score, ind_t] = nanmin(score,[],2);
            % m_dist = pdistCS(ind_t);
            
            % Then convert m_socre to phi_classifier (phi_score)
            C = 7;  % C is slope
            H = 1;  % H is the truncate position for mDist
            % With this, phi_classifier = (SF - 7*mDist)/8 + 15/16
            phi_score = stru(iS).tSF(ind_t)/(1+C*H) -C/(1+C*H)*min(H,m_dist) + (1+2*C*H)/2/(1+C*H);
            
            tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
            if ~isempty(tNote.enable) && ismember([ID_current,iCluster],tNote.enable,'rows')
                stru(iS).cEnable(iCluster) = 1;
                
            end
            if ~isempty(tNote.disable) && ismember([ID_current,iCluster],tNote.disable,'rows')
                if (phi_score > scoreCF)
                    stru(iS).cEnable(iCluster) = -1;        % to qualify to 'disable', it has to be identified as twin if were not disabled (phi_score < scoreCF)
                end
            end
            
            
            if tsNum > nss
                sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
                mDistMapLocal(indClusterLocal) = m_dist;
                scoreMapLocal(indClusterLocal) = phi_score;

                if ((phi_score > scoreCF)&&(-1 ~= stru(iS).cEnable(iCluster))) || (1 == stru(iS).cEnable(iCluster))
                    if phi_score > scoreCF
                        % (1) condition: if phi_score < scoreCF, then no need to use 'Enable', so set the .cEnable(iCluster) back to 0
                        stru(iS).cEnable(iCluster) = 0;
                        scoreMapLocal(indClusterLocal) = phi_score + 1;
                    else
                        % (2) condition: else, it needs to be enabled
                        scoreMapLocal(indClusterLocal) = scoreForEnabled;
                    end
                    twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map. For twinMap, assign if phi_score < scoreCF
                    stru(iS).c2t(iCluster) = tsNum;     % c2t is the identification label.  Cluster->Twin
                end
                if (-1 == stru(iS).cEnable(iCluster))
                    % (3) condition: if disabled
                    twinMapLocal(indClusterLocal) = -tsNum;
                    scoreMapLocal(indClusterLocal) = scoreForDisabled;
                    stru(iS).c2t(iCluster) = tsNum;
                end
                
            end
            
        end
        
        % copy identified twin system number to twinMap
        twinMap(indR_min:indR_max, indC_min:indC_max) = twinMapLocal;
        sfMap(indR_min:indR_max, indC_min:indC_max) = sfMapLocal;
        mDistMap(indR_min:indR_max, indC_min:indC_max) = mDistMapLocal;
        scoreMap(indR_min:indR_max, indC_min:indC_max) = scoreMapLocal;

    end
    waitbar(iS/length(stru), hWaitbar);
    %     input('press to continue');
end

scoreMap(scoreMap==0)=nan;

try
    close(hWaitbar);
catch
end

if saveTF
    timeStr = datestr(now,'yyyymmdd_HHMM');
    save(['temp_tNote_s',num2str(iE),'_',timeStr,'.mat'],'tNote','scoreCF');
    save([saveDataPath,'twin_label_tNote_s',num2str(iE),'_',timeStr,'.mat'],'iE','tNote','scoreCF');
end

%% save the result (including stru, tNote, maps) at the end
name_cluster_to_twin_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];

save([saveDataPath,name_cluster_to_twin_result],'clusterNumMap','twinMap','sfMap','mDistMap','scoreMap','stru','tNote','scoreCF','-append');
try
    save([saveDataPath,name_cluster_to_twin_result],'clusterNumMapCleaned','-append');
end
% save([saveDataPath,name_result_modified],'stru','-append');
% save([saveDataPath,name_result_modified],'twinMap_2','shearMap','sfMap_2','costMap','-append');
disp('saved cluster_to_twin_result');


% temp code for plot and investigate the results
if 0
    %     myplot(X,Y,exx,boundaryTFB);
    %     myplot(X,Y,clusterNumMap,boundaryTFB);
    %
    %     myplot(X,Y,mDistMap,boundaryTFB);
    %     myplot(X,Y,twinMap,boundaryTFB); caxis([18,24]);
    %     myplot(X,Y,sfMap,boundaryTFB);
    %
    %     myplot(X,Y,shearMap, boundaryTFB);
    %     myplot(X,Y,abs(shearMap-0.1289), boundaryTFB); title('shear diff');
    %     myplot(X,Y,twinMap_2, boundaryTFB); caxis([18,24]);
    %     myplot(X,Y,sfMap_2, boundaryTFB);
    %     myplot(X,Y,costMap, boundaryTFB);
end


%% This loads the strain maps at all strain levels considered.
for iE = iE_start:iE_stop
    name_cluster_to_twin_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,name_cluster_to_twin_result],'clusterNumMap');
    cnMap{iE} = clusterNumMap;
    
    strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile)
    clear('exy_corrected');
    load(strainFile,'exx','exy_corrected');     % if 'exy_corrected' does not exist, this does not give error, rather, just warning. % ----------------------------------------------------------------------------------
    %     load(strainFile,'exx','exy','eyy','sigma','exy_corrected');
    if exist('exy_corrected','var')&&(1==exy_corrected)
        disp('================= exy already corrected ! ========================');
        exy_corrected = 1;
    else
        disp('================= exy being corrected here ! =======================');
        exy = -exy;
        exy_corrected = 1;
    end
    % remove bad data points
    exx(sigma==-1) = nan;
    %     exy(sigma==-1) = nan;
    %     eyy(sigma==-1) = nan;
    qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
    %     qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
    %     qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
    ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
    exx(ind_outlier) = nan;
    %     exy(ind_outlier) = nan;
    %     eyy(ind_outlier) = nan;
    es{iE} = exx;
end
