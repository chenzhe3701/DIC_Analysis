% Start with analyzing more data, such as quantiles, ...

clear;
addChenFunction;
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);

gIDwithTrace = gID(~isnan(gExx));
% gIDwithTrace = [89,129,135,191,201,210,327,401,422,553];        % WE43 T5 #7, new
% gIDwithTrace = [302,151,186,191,296,431,572,1211];            % WE43_T6_C1

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 6;
resReduceRatio = 3;         % to save space, reduce map resolution
grow_boundary_TF = 0;       % whether to grow boundary to make it thicker
% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'


%% Can load strain data for a specific strain level
rng(1);
for iE = iE_start:iE_stop
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
    
    struNC_fName = [sampleName, '_s',STOP{iE+B},'_nClusters.mat']
    load([saveDataPath,struNC_fName],'struNC');
    %% Cluster strain and record information in a structure, assign clusters to twin systems on the fly, but can also modify later.
    
    debugTF = 0;
    %     distCI = 0.035;     % distance Criterion Initial
    sfCI = -0.5;        % schmid factor Criterion Initial
    shearTarget = 0.1289;
    shearCI = 0.05;        % shear criterion initial
    costCI = 0.07;     % sqrt(cost) criterion initial
    scoreCI = 0.19;          % score criterion initial, score = 2*sqrt(cost) + abs(shearFit-shearTarget)
    
    % Create a few maps to record the criterion.
    clusterNumMap = zeros(size(exx));
    twinMap = zeros(size(exx));
    sfMap = zeros(size(exx));
    disSimiMap = zeros(size(exx));
    scoreMap = zeros(size(exx));
    
    twinMap_2 = zeros(size(exx));
    shearMap = zeros(size(exx));
    sfMap_2 = zeros(size(exx));
    costMap = zeros(size(exx));
    
    [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
    ss = crystal_to_cart_ss(ssa,c_a);
    
    % gIDwithTrace = [86, 41, 193, 194, 259, 262, 197, 296, 153, 182, 378, 451, 576]; % select some obviously twinned grain for debugging
    
    %%
    clear stru;
    hWaitbar = waitbar(0,'finding twin region for grains ...');
    for iS =1:length(gIDwithTrace)
        %%
        %         iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        %         iS = find(gIDwithTrace == 296); % for debugging.
        try
            stru(iS)=[];
        catch
        end
        close all;
        % select the target grain
        ID_current = gIDwithTrace(iS);  % id=262 for an example for WE43-T6-C1
        
        % =================== find measured strain: =======================
        ind_current = find(ID_current == gID);    % an index of row
        ID_neighbor = gNeighbors(ind_current,:);
        ID_neighbor = ID_neighbor((ID_neighbor~=0)&(ID_neighbor~=neighbor_elim));
        
        % find index range of a small matrix containing the grain of interest,
        % (and can choose to include some neighboring grains)
        ind_local = ismember(ID, [ID_current]); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
        exy_local = exy(indR_min:indR_max, indC_min:indC_max);
        eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        x_local = X(indR_min:indR_max, indC_min:indC_max);
        y_local = Y(indR_min:indR_max, indC_min:indC_max);
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        
        % find vectors for cluster, using ind
        ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));
        exx_t = exx_local(ind);
        exy_t = exy_local(ind);
        eyy_t = eyy_local(ind);
        data_t = [exx_t(:), exy_t(:), eyy_t(:)];
        %     [data_zn,mean0,std0] = zero_normalize_column(data_t);   % seems like it's better without zero_normalize.
        
        % ============= clustering data.  grain-197 is a good example showing that kmeans seems to be better than gmModels =====================
        % first, predict centroid
        stru(iS).gID = ID_current;
        
        % ================ method-1, from theoretical twin shear --> to predict strain components as cluster center ===================
        ind_euler = find(gID==ID_current);
        euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
        if (1==eulerAligned)
            g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
        else
            g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
        end
        gamma = 0.1289; % twin shear for Mg
        cPred = nan*zeros(nss,5);   % [iss, SF, exx, exy, eyy]
        for iss = (nss+1):(nss+ntwin)   % for Mg
            %         disp('---');
            N(iss,:) = ss(1,:,iss) * g;
            M(iss,:) = ss(2,:,iss) * g;
            MN2{iss} = M(iss,:)'*N(iss,:);
            MN2{iss} = MN2{iss}(1:2,1:2);
            %         F3 = eye(3) + gamma*M(iss,:)'*N(iss,:);
            %         F = F3(1:2,1:2);
            F = eye(2) + gamma*MN2{iss};
            epsilon = (F'*F-eye(2))/2;
            %         disp((F3'*F3-eye(3))/2);
            %         disp(epsilon);
            cPred(iss,1) = iss;                                     % ss number
            cPred(iss,2) = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
            cPred(iss,3:5) = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.  Note that 'conjugated' twin system, i.e., 19 and 22, almost always show similar components!!!
        end
        stru(iS).tLabel = (nss+1 : nss+ntwin)';         % twin system number
        stru(iS).tSF = cPred(nss+1:nss+ntwin,2)';       % twin schmid factor
        stru(iS).tStrain = cPred(nss+1:nss+ntwin,3:5);      % twin strain components
        [~,ind_centroid_initial] = max(stru(iS).tSF);
        centroid_initial = stru(iS).tStrain(ind_centroid_initial,:);
        %     disp(cPred);
        
        if(~isempty(struNC(iS).nCluster))
            ii = find(arrayfun(@(x) x.gID == stru(iS).gID, struNC));
            nCluster = struNC(ii).nCluster;
            
            % ========================= (1) perform kmeans cluster ============================
            %         nCluster = 4;       % total number of clusters
            nRep = 3;
            c0 = kmeans_pp_init(data_t,nCluster,nRep,centroid_initial);
            [idx, centroid, sumd] = kmeans(data_t, nCluster, 'Distance','sqeuclidean','MaxIter',1000,'start',c0);   % 'correlation' distance not good.
            clusterNumMapLocal = zeros(size(x_local));      % record raw clusterNumberMap
            clusterNumMapLocal(ind) = idx;
            clusterNumMap(indR_min:indR_max, indC_min:indC_max) = clusterNumMap(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
            
            stru(iS).cLabel = (1:nCluster)';     % cluster number (actually label, but labe=number)
            stru(iS).cCen = centroid;          % cluster centroid
            
            % ==================== continue method-1, assign cluster to twin system on the fly ==========================
            twinMapLocal = zeros(size(exx_local));          % local map to record twin_system_number
            sfMapLocal = zeros(size(exx_local));            % local map to record schmid_factor
            disSimiMapLocal = zeros(size(exx_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
            scoreMapLocal = zeros(size(exx_local));
            
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                indClusterLocal = (clusterNumMapLocal==cNum);
                
                % [For debug], extract local exx,exy,eyy map of this cluster, and local data [exx(:),exy(:),eyy(:)], for whatever analysis -----------------  
                exx_cluster = zeros(size(exx_local));   % map
                exx_cluster(indClusterLocal) = exx_local(indClusterLocal);
                exy_cluster = zeros(size(exy_local));
                exy_cluster(indClusterLocal) = exy_local(indClusterLocal);
                eyy_cluster = zeros(size(eyy_local));
                eyy_cluster(indClusterLocal) = eyy_local(indClusterLocal);
                data_cluster = [exx_cluster(indClusterLocal),exy_cluster(indClusterLocal),eyy_cluster(indClusterLocal)];    % matrix
                
                egcName = [num2str(iE*100000+ID_current*10+iCluster)];  % name = strain-grain-cluster, for saving purpose
                
                % [For debug end] whatever------------------------------------------------------------------------------------  
                
                
                pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
                %             pdistCS(pdistCS > distCI) = nan;                 % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
                %             pdistCS(stru(iS).tSF < sfCI) = nan;             % [criterion-2] SF must > sfCI
                %             [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
                
                % score boundary y=kx+b passes [pdist2,sf] = [0, 0.15] and [0.05, 0.5]
                score = pdistCS * 7 + (0.5 - stru(iS).tSF);
                [m_score, ind_t] = nanmin(score,[],2);
                m_dist = pdistCS(ind_t);
                
                tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
                
                if tsNum > -1
                    % Note global assign like: twinMap(ind_c_on_map) = tsNum; is very slow.
                    twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map
                    sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
                    disSimiMapLocal(indClusterLocal) = m_dist;
                    scoreMapLocal(indClusterLocal) = m_score;
                end
            end
            if 1 == debugTF
                
                myplot(x_local,y_local,exx_local); % myplot(x_local,y_local,exy_local); myplot(x_local,y_local,eyy_local);
                myplot(x_local,y_local,clusterNumMapLocal);
                % myplot(x_local,y_local,twinMapLocal); caxis([nss,nss+ntwin]);
                % myplot(x_local,y_local,sfMapLocal);
                % myplot(x_local,y_local,disSimiMapLocal);
            end
            % copy identified twin system number to twinMap
            twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
            sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
            disSimiMap(indR_min:indR_max, indC_min:indC_max) = disSimiMap(indR_min:indR_max, indC_min:indC_max) + disSimiMapLocal;
            scoreMap(indR_min:indR_max, indC_min:indC_max) = scoreMap(indR_min:indR_max, indC_min:indC_max) + scoreMapLocal;
            
            
            
            % ============ method-2, from cluster_centroid to estimated_shear, and compare if this is similar enough to theoretical shear of 0.1289 ====================
            options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');
            clear shear cost;
            for iCluster=1:nCluster
                for iss = (nss+1):(nss+ntwin)
                    [shear(iCluster,iss),cost(iCluster,iss)]=fminunc(@(x) sqrt(sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
                        [centroid(iCluster,1),centroid(iCluster,2);centroid(iCluster,2),centroid(iCluster,3)]).^2))), 0.1298, options);
                end
            end
            stru(iS).cShear = shear(:,nss+1:nss+ntwin);       % cluster centroid's expected shear
            stru(iS).cCost = cost(:,nss+1:nss+ntwin);      % cluster centroid fitted shear -induced cost
            
            % =============== continue method-2, assign cluster to twin system on the fly ==============================
            twinMapLocal_2 = zeros(size(exx_local));          % local map to record twin_system_number
            shearMapLocal = zeros(size(exx_local));            % local map to record schmid_factor
            sfMapLocal_2 = zeros(size(exx_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
            costMapLocal = zeros(size(exx_local));       % local map to record cost
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                indClusterLocal = (clusterNumMapLocal==cNum);
                
                shearFit = stru(iS).cShear(iCluster,:);
                shearFit(stru(iS).tSF < sfCI) = nan;            % [criterion-1] SF must > 0.35
                shearFit(abs(shearFit-shearTarget) > shearCI) = nan;      % [criterion-2] shear_difference should be < 0.5
                
                costFit = stru(iS).cCost(iCluster,:);
                costFit(costFit > costCI) = nan;     % [criterion-3] costFit should be < 0.025.  Change to 0.035 empirically
                
                score = abs(shearFit-shearTarget)+2*costFit;     % [criterion-4] the criterion is a combination of shearFit and sqrt(cost)
                score(costFit*2+abs(shearFit-shearTarget) > scoreCI) = nan;  % [criterion-5] combination of shearDiff and sqrtCost should also be limited
                
                [m_score,ind_t] = min(score,[],2);
                % [m_shear_diff,ind_t] = min(abs(shearFit-0.1289),[],2);    % if just look at the shear.
                
                if isnan(m_score)
                    tsNum = [];
                else
                    tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
                end
                if tsNum > nss
                    twinMapLocal_2(indClusterLocal) = tsNum;
                    sfMapLocal_2(indClusterLocal) = stru(iS).tSF(ind_t);
                    shearMapLocal(indClusterLocal)  = stru(iS).cShear(iCluster,ind_t);
                    costMapLocal(indClusterLocal) = stru(iS).cCost(iCluster,ind_t);
                end
            end
            if 1 == debugTF
                myplot(x_local,y_local,twinMapLocal_2); caxis([nss,nss+ntwin]);
                myplot(x_local,y_local,sfMapLocal_2);
                myplot(x_local,y_local,shearMapLocal);
                myplot(x_local,y_local,costMapLocal);
            end
            % copy identified twin system number to twinMap
            twinMap_2(indR_min:indR_max, indC_min:indC_max) = twinMap_2(indR_min:indR_max, indC_min:indC_max) + twinMapLocal_2;
            shearMap(indR_min:indR_max, indC_min:indC_max) = shearMap(indR_min:indR_max, indC_min:indC_max) + shearMapLocal;
            sfMap_2(indR_min:indR_max, indC_min:indC_max) = sfMap_2(indR_min:indR_max, indC_min:indC_max) + sfMapLocal_2;
            costMap(indR_min:indR_max, indC_min:indC_max) = costMap(indR_min:indR_max, indC_min:indC_max) + costMapLocal;
            
            try
                waitbar(iS/length(gIDwithTrace), hWaitbar);
            catch
            end
            
            if 1==debugTF
                temp=stru(iS)
                disp(['iS=',num2str(iS)]);
                %         input('press to continue');
            end
            disp(['finished ID=',num2str(ID_current)]);
            
        end
    end
    
    scoreMap(scoreMap==0)=nan;
    
    try
        close(hWaitbar);
    catch
    end
    %%
    timeStr = datestr(now,'yyyymmdd_HHMM')
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly_',timeStr,'.mat'];
    save([saveDataPath,name_result_on_the_fly] ,'clusterNumMap','twinMap','sfMap','disSimiMap','scoreMap','twinMap_2','shearMap','sfMap_2','costMap','stru');
    disp(['saved cluster_result_on_the_fly: ',name_result_on_the_fly]);
    
    
end





