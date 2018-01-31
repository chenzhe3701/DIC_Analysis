% Start with analyzing more data, such as quantiles, ...

clear;
addChenFunction;
dicPath = uigetdir('','pick DIC directory, which contains the stitched DIC data for each stop');
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
iE_start = 1;   % elongation levels to analyze. 0-based.
iE_stop = 7;
resReduceRatio = 3;         % to save space, reduce map resolution
grow_boundary_TF = 0;       % whether to grow boundary to make it thicker
% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

notes = struct('atEdge',[],'likeTwin',[],'tooSmall',[]);
% This is a special note for data processing, for specific samples, to give grain labels
% useful fields will be 'tooSmall', 'atEdge', 'likeTwin'
% notes = load('T5#7_traceAnalysisNotes');
% 2016-08-15. create a field called 'noTrace'. noTrace = 1 = too small. noTrace = 2 = at edge .

% end of modify settings part 1 ------------------------------------------------------------------------------------------------------------------------------------
save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],...
    'dicPath','dicFiles',...
    'STOP','iE_start','iE_stop','resReduceRatio','grow_boundary_TF','f1','f2','neighbor_elim','twinTF_text','notes','gIDwithTrace',...
    '-append');

%%  can load stru
% iE = 5;
% name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
% try
%     load([saveDataPath,name_result_on_the_fly]);
% catch
% end
%% Can load strain data for a specific strain level
for iE = iE_start:iE_stop
    strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
    load(strainFile,'exx','exy','eyy');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------
    
    
    %% Cluster strain and record information in a structure, assign clusters to twin systems on the fly, but can also modify later.
    
    debugTF = 0;
    distCI = 0.05;     % distance Criterion Initial
    sfCI = 0.3;        % schmid factor Criterion Initial
    shearTarget = 0.1289;
    shearCI = 0.05;        % shear criterion initial
    sqrtCostCI = 0.06;     % sqrt(cost) criterion initial
    scoreCI = 0.13;          % score criterion initial, score = 2*sqrt(cost) + abs(shearFit-shearTarget)
    
    % Create a few maps to record the criterion.
    clusterNumMap = zeros(size(exx));
    twinMap = zeros(size(exx));
    sfMap = zeros(size(exx));
    disSimiMap = zeros(size(exx));
    twinMap_2 = zeros(size(exx));
    shearMap = zeros(size(exx));
    sfMap_2 = zeros(size(exx));
    costMap = zeros(size(exx));
    
    [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
    ss = crystal_to_cart_ss(ssa,c_a);
    
    % gIDwithTrace = [86, 41, 193, 194, 259, 262, 197, 296, 153, 182, 378, 451, 576]; % select some obviously twinned grain for debugging
    % tune up grains, step 5:
    % 129, 354, 451, 577, 691, 783, 797, 1141, 1266, 1280, 1350,
    % 110, 195, 569, 1266, 1313, 1141
    % step 2: 195, 541
    
    %%
    hWaitbar = waitbar(0,'finding twin region for grains ...');
    for iS =1:length(gIDwithTrace)
        %         iS = find(arrayfun(@(x) x.gID == 1466,stru)),  % for debugging
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
        
        nRow = indR_max - indR_min + 1;
        nColumn = indC_max - indC_min + 1;
        
        exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
        exy_local = exy(indR_min:indR_max, indC_min:indC_max);
        eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        x_local = X(indR_min:indR_max, indC_min:indC_max);
        y_local = Y(indR_min:indR_max, indC_min:indC_max);
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        % (0) can smooth
        exx_local = colfilt(exx_local, [3 3], 'sliding', @(x) nanmean(x,1));
        exy_local = colfilt(exy_local, [3 3], 'sliding', @(x) nanmean(x,1));
        eyy_local = colfilt(eyy_local, [3 3], 'sliding', @(x) nanmean(x,1));
        
        % find vectors for cluster, using ind
        ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));
        exx_t = exx_local(ind);
        exy_t = exy_local(ind);
        eyy_t = eyy_local(ind);
        data_t = [exx_t(:), exy_t(:), eyy_t(:)];
        %     [data_zn,mean0,std0] = zero_normalize_column(data_t);   % seems like it's better without zero_normalize.
        
        % export data if needed
        %         save(['grain_',num2str(ID_current)],'data_t');
        
        % determine optimum number of clusters
        maxCluster = 6;
        nPoints = 2000;
        ind_reduce = ~isnan(sum(data_t,2));
        data_reduce = data_t(ind_reduce,:);
        reduce_ratio = fix(size(data_reduce,1)/nPoints);
        data_reduce = data_reduce(1:reduce_ratio:end,:);
        
        if isempty(data_reduce)
            nCluster = 4;
        else
            eva = evalclusters(data_reduce,'kmeans','silhouette','klist',2:maxCluster)
            nCluster = eva.OptimalK;
        end
        
        % compare, test the 'stability' of silhouette: calculate from the result of another kmeans run, and compare.
        %         clear withinSum score;
        %         for nc = 2:maxCluster
        %             [idx_cell{nc}, centroid, sumd_cell{nc}] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000);   % 'correlation' distance not good.
        %             [score_cell{nc}] = silhouette(data_reduce,idx_cell{nc});
        %             withinSum(nc) = mean(sumd_cell{nc});
        %             score(nc) = nanmean(score_cell{nc});
        %             figure; silhouette(data_reduce,idx_cell{nc});
        %         end
        
        %         figure;
        %         subplot(1,2,1);
        %         plot(2:maxCluster,withinSum(2:end),'x-'); xlabel('num of clusters'); ylabel('within ssd'); axis square; % within-cluster SS, always decrease, maybe no need to look at.
        
        %         figure;
        %         try
        %             subplot(1,2,1);
        %             plot(eva);
        %             axis square;
        %             subplot(1,2,2);
        %         catch
        %             subplot(1,1,1);
        %         end
        %         plot(2:maxCluster,score(2:end),'x-'); xlabel('num of clusters'); ylabel('avg silhouette'); axis square;
        
        
        % compare 4 criterions
        %         for nPoints = [1000,2000,5000,10000]
        %             ind_reduce = ~isnan(sum(data_t,2));
        %             data_reduce = data_t(ind_reduce,:);
        %             reduce_ratio = fix(size(data_reduce,1)/nPoints);
        %             data_reduce = data_reduce(1:reduce_ratio:end,:);
        %             eva1 = evalclusters(data_reduce,'kmeans','silhouette','klist',2:maxCluster)
        %             figure; subplot(2,2,1); plot(eva1);
        %             eva2 = evalclusters(data_reduce,'kmeans','gap','klist',2:maxCluster)
        %             subplot(2,2,2); plot(eva2);
        %             eva3 = evalclusters(data_reduce,'kmeans','CalinskiHarabasz','klist',2:maxCluster)
        %             subplot(2,2,3); plot(eva3);
        %             eva4 = evalclusters(data_reduce,'kmeans','DaviesBouldin','klist',2:maxCluster)
        %             subplot(2,2,4); plot(eva4);
        %         end
        
        
        %% ============= clustering data.  grain-197 is a good example showing that kmeans seems to be better than gmModels =====================
        % (1) kmeans cluster
        nCluster = 2;       % total number of clusters
        [idx, centroid, sumd] = kmeans(data_t, nCluster, 'Distance','sqeuclidean','MaxIter',1000,'replicates',1);   % 'correlation' distance not good.
        clusterNumMapLocal = zeros(size(x_local));      % record raw clusterNumberMap
        clusterNumMapLocal(ind) = idx;
        clusterNumMap(indR_min:indR_max, indC_min:indC_max) = clusterNumMap(indR_min:indR_max, indC_min:indC_max) + clusterNumMapLocal;
        
        
        %     % try to solve equation, hasn't been able to
        %     syms shear
        %     vpasolve( 1/2*((eye(2)+shear*MN{19}(1:2,1:2))'*(eye(2)+shear*MN{19}(1:2,1:2))-eye(2)) == [centroid(1,1), centroid(1,2); centroid(1,2), centroid(1,3)],shear)
        
        %     % (1.1) zero-normalize, then kmeans cluster
        %     [idx, centroid, sumd] = kmeans(data_zn,4,'Distance','sqeuclidean','MaxIter',1000);   % 'correlation' distance not good.
        %     clusterNum = zeros(size(x_local,1),size(x_local,2));
        %     clusterNum(ind) = idx;
        %     fh11 = myplot(x_local,y_local,clusterNum);
        
        
        %     % (2) hierarchical cluster.  Does not work well.
        %     Z = linkage(data_t);
        %     figure;dendrogram(Z);
        %     T = cluster(Z,'maxclust',3);
        %     clusterNum = zeros(size(x_local,1),size(x_local,2));
        %     clusterNum(ind) = T;
        %     fh2 = myplot(x_local,y_local,clusterNum);
        
        %     % (3) using GMM Gaussian Mixture Models, not as good as kmeans
        %     AIC = zeros(1,4);
        %     gmModels = cell(1,4);
        %     options = statset('MaxIter',1000);
        %     for k = 1:4
        %         gmModels{k} = fitgmdist(data_t,k,'options',options,'SharedCovariance',false);
        %         AIC(k) = gmModels{k}.AIC;
        %     end
        %     [minAic,numComponents] = min(AIC);
        %
        %     T = cluster(gmModels{numComponents}, data_t);
        %     clusterNum = zeros(size(x_local,1),size(x_local,2));
        %     clusterNum(ind) = T;
        %     fh3 = myplot(x_local,y_local,clusterNum);
        
        stru(iS).gID = ID_current;
        stru(iS).cLabel = (1:nCluster)';     % cluster number (actually label, but labe=number)
        stru(iS).cCen = centroid;          % cluster centroid
        
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
        %     disp(cPred);
        
        % ==================== continue method-1, assign cluster to twin system on the fly ==========================
        twinMapLocal = zeros(size(exx_local));          % local map to record twin_system_number
        sfMapLocal = zeros(size(exx_local));            % local map to record schmid_factor
        disSimiMapLocal = zeros(size(exx_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
        
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            indClusterLocal = (clusterNumMapLocal==cNum);
            %  ind_c_on_map = (ID==ID_current)&(clusterNumMap == cNum);  % could index on whole map, but very slow
            
            pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
            pdistCS(pdistCS > distCI) = nan;                 % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
            pdistCS(stru(iS).tSF < sfCI) = nan;             % [criterion-2] SF must > sfCI
            [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
            
            if isnan(m_dist)
                tsNum = [];
            else
                tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
            end
            if tsNum > nss
                % Note global assign like: twinMap(ind_c_on_map) = tsNum; is very slow.
                twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map
                sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
                disSimiMapLocal(indClusterLocal) = m_dist;
            end
        end
        if 1 == debugTF
            myplot(x_local,y_local,exx_local);
            %          myplot(x_local,y_local,exy_local);
            %          myplot(x_local,y_local,eyy_local);
            myplot(x_local,y_local,clusterNumMapLocal);
            %             myplot(x_local,y_local,twinMapLocal); caxis([nss,nss+ntwin]);
            %             myplot(x_local,y_local,sfMapLocal);
            %             myplot(x_local,y_local,disSimiMapLocal);
        end
        % copy identified twin system number to twinMap
        twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
        sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
        disSimiMap(indR_min:indR_max, indC_min:indC_max) = disSimiMap(indR_min:indR_max, indC_min:indC_max) + disSimiMapLocal;
        
        % ============ method-2, from cluster_centroid to estimated_shear, and compare if this is similar enough to theoretical shear of 0.1289 ====================
        options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');
        for iCluster=1:nCluster
            for iss = (nss+1):(nss+ntwin)
                [shear(iCluster,iss),cost(iCluster,iss)]=fminunc(@(x) sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
                    [centroid(iCluster,1),centroid(iCluster,2);centroid(iCluster,2),centroid(iCluster,3)]).^2)), 0.1298, options);
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
            
            shearFit = stru(iS).cShear(iCluster,:);
            shearFit(stru(iS).tSF < sfCI) = nan;            % [criterion-1] SF must > 0.35
            shearFit(abs(shearFit-shearTarget) > shearCI) = nan;      % [criterion-2] shear_difference should be < 0.5
            
            costFit = stru(iS).cCost(iCluster,:);
            costFit(sqrt(costFit) > sqrtCostCI) = nan;     % [criterion-3] costFit should be < 0.025.  Change to 0.035 empirically
            
            score = abs(shearFit-shearTarget)+2*sqrt(costFit);     % [criterion-4] the criterion is a combination of shearFit and sqrt(cost)
            score(sqrt(costFit)*2+abs(shearFit-shearTarget) > scoreCI) = nan;  % [criterion-5] combination of shearDiff and sqrtCost should also be limited
            
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
    end
    close(hWaitbar);
    %%
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    disp('start saving result_on_the_fly');
    save([saveDataPath,name_result_on_the_fly] ,'clusterNumMap','twinMap','sfMap','disSimiMap','twinMap_2','shearMap','sfMap_2','costMap','stru');
    disp('finished saving result_on_the_fly');
    
    
    %% =========== match cluster with twin system again, if need to change parameter ============================
    iE = iE;
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    load([saveDataPath,name_result_on_the_fly],'clusterNumMap');
    
    % criterion final, for modifying
    distCF = 0.035;
    sfCF = 0.37;
    shearTarget = 0.1289;
    shearCF = 0.05;
    sqrtCostCF = 0.042;
    scoreCF = 0.1;
    
    close all;
    % Create a few maps to record the criterion.
    twinMap = zeros(size(exx));
    sfMap = zeros(size(exx));
    disSimiMap = zeros(size(exx));
    twinMap_2 = zeros(size(exx));
    shearMap = zeros(size(exx));
    sfMap_2 = zeros(size(exx));
    costMap = zeros(size(exx));
    
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
        
        twinMapLocal_2 = zeros(size(ID_local));          % local map to record twin_system_number
        shearMapLocal = zeros(size(ID_local));            % local map to record schmid_factor
        sfMapLocal_2 = zeros(size(ID_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
        costMapLocal = zeros(size(ID_local));       % local map to record cost
        
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            indClusterLocal = (clusterNumMapLocal==cNum);
            
            % ============== method-1 =============================
            pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
            pdistCS(pdistCS > distCF) = nan;                 % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
            pdistCS(stru(iS).tSF < sfCF) = nan;             % [criterion-2] SF must > 0.35
            [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
            
            if isnan(m_dist)
                tsNum = [];
            else
                tsNum = stru(iS).tLabel(ind_t);               % match cluster to this twin system
            end
            if tsNum > nss
                twinMapLocal(indClusterLocal) = tsNum;    % assign twinSysNum to the region in the local map
                sfMapLocal(indClusterLocal) = stru(iS).tSF(ind_t);
                disSimiMapLocal(indClusterLocal) = m_dist;
            end
            
            % ================ method-2 ==========================
            shearFit = stru(iS).cShear(iCluster,:);
            shearFit(stru(iS).tSF < sfCF) = nan;            % [criterion-1] SF must > 0.35
            shearFit(abs(shearFit-shearTarget) > shearCF) = nan;      % [criterion-2] shear_difference should be < 0.5
            
            costFit = stru(iS).cCost(iCluster,:);
            costFit(sqrt(costFit) > sqrtCostCF) = nan;                 % [criterion-3] costFit should be < 0.025.  Change to 0.035 empirically
            
            score = abs(shearFit-shearTarget) + 2*sqrt(costFit);       % [criterion-4] the criterion is a combination of shearFit and sqrt(cost)
            score(sqrt(costFit)*2+abs(shearFit-shearTarget) > scoreCF) = nan;      % [criterion-5] combination of shearDiff and sqrtCost should also be limited
            
            [m_score,ind_t] = min(score,[],2);
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
        
        % copy identified twin system number to twinMap
        twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
        sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
        disSimiMap(indR_min:indR_max, indC_min:indC_max) = disSimiMap(indR_min:indR_max, indC_min:indC_max) + disSimiMapLocal;
        
        twinMap_2(indR_min:indR_max, indC_min:indC_max) = twinMap_2(indR_min:indR_max, indC_min:indC_max) + twinMapLocal_2;
        shearMap(indR_min:indR_max, indC_min:indC_max) = shearMap(indR_min:indR_max, indC_min:indC_max) + shearMapLocal;
        sfMap_2(indR_min:indR_max, indC_min:indC_max) = sfMap_2(indR_min:indR_max, indC_min:indC_max) + sfMapLocal_2;
        costMap(indR_min:indR_max, indC_min:indC_max) = costMap(indR_min:indR_max, indC_min:indC_max) + costMapLocal;
        
        waitbar(iS/length(stru), hWaitbar);
        %     input('press to continue');
    end
    close(hWaitbar);
    
    name_result_modified = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_modified.mat'];
    disp('start saving cluster_result_modified');
    save([saveDataPath,name_result_modified],'twinMap','sfMap','disSimiMap','twinMap_2','shearMap','sfMap_2','costMap','clusterNumMap','stru');
    disp('finished saving cluster_result_modified');
    
    
    %% temp code for plot and investigate the results
    if 0
        myplot(X,Y,exx,boundaryTFB);
        myplot(X,Y,clusterNumMap,boundaryTFB);
        
        myplot(X,Y,disSimiMap,boundaryTFB);
        myplot(X,Y,twinMap,boundaryTFB); caxis([18,24]);
        % myplot(X,Y,sfMap,boundaryTFB);
        
        myplot(X,Y,shearMap, boundaryTFB);
        myplot(X,Y,abs(shearMap-0.1289), boundaryTFB); title('shear diff');
        myplot(X,Y,twinMap_2, boundaryTFB); caxis([18,24]);
        myplot(X,Y,sfMap_2, boundaryTFB);
        myplot(X,Y,costMap, boundaryTFB);
    end
    
end
%%
% for iS = 1%:length(gIDwithTrace)
%     close all;
%     ID_current = gIDwithTrace(iS);              % id of current grain
%     ID_current = 668;
%
%     ind_current = find(ID_current == gID);    % an index of row
%     phi1_current = gPhi1(ind_current);
%     phi_current = gPhi(ind_current);
%     phi2_current = gPhi2(ind_current);
%
%     ID_neighbor = gNeighbors(ind_current,:);
%     ID_neighbor = ID_neighbor((ID_neighbor~=0)&(ID_neighbor~=neighbor_elim));
%
%     % find index range of a small matrix containing the grain of interest
%     ind_pool = ismember(ID, [ID_current,ID_neighbor]);
%     indC_min = find(sum(ind_pool, 1), 1, 'first');
%     indC_max = find(sum(ind_pool, 1), 1, 'last');
%     indR_min = find(sum(ind_pool, 2), 1, 'first');
%     indR_max = find(sum(ind_pool, 2), 1, 'last');
%
%     nRow = indR_max - indR_min + 1;
%     nColumn = indC_max - indC_min + 1;
%
%     e_current = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
%     boundaryTF_current = boundaryTF(indR_min:indR_max, indC_min:indC_max);
%     x_current = X(indR_min:indR_max, indC_min:indC_max);
%     y_current = Y(indR_min:indR_max, indC_min:indC_max);
%     ID_map_current = ID(indR_min:indR_max, indC_min:indC_max);
%
%     e_grain = e_current;
%     e_grain(ID_map_current~=ID_current) = 0;  % 'e_grain' is strain of This grain. 'e_current' is strian of this region.
%     e_grain(isnan(e_grain)) = 0;
%
%     % calculate Schmid factor.
%     % sf_mat = [#, SF, angle_XtoY, trace_x_end, trace_y_end]
%     % burgersXY = [burgers_X, burgers_Y, ratio].
%     % ---------------------------------------- select proper setting for analysis, such as material, twin, stress ---------------------------------------------------------
%
%     [sf_mat, sf_mat_sorted, burgersXY] = trace_analysis_TiMgAl([phi1_current,phi_current,phi2_current],[-90,180,0],[0,0,0],stressTensor,sampleMaterial,twinTF_text);
%
%
%
%     figureHandle_1 = myplot(x_current,y_current,e_grain);
%
%     % hough transform, houghpeaks the peaks, houghlines the line segments
%     [h,t,r] = hough(e_grain);
%     [tg,rg] = meshgrid(t,r);
%     % figureHandle_2 = myplot(tg,rg,h);
%     figureHandle_2 = figure;
%     imshow(h,[],'XData',t,'YData',r,'colormap',parula); axis on; axis square;
%
%     peaks = houghpeaks(h,10);   %,'NHoodSize',[floor(size(h,1)/2)*2+1,3]);
%     set(0,'currentfigure',figureHandle_2); hold on;
%     for k = 1:size(peaks,1)
%         xy = [tg(peaks(k,1),peaks(k,2)),rg(peaks(k,1),peaks(k,2))];
%         plot3(xy(1),xy(2),max(h(:)),'s','LineWidth',1,'Color','k')
%     end
%
%     lines = houghlines(e_grain,t,r,peaks);
%     % lines(k).point1/2 is in fact the index (index_c, index_r)
%     % show the extracted lines
%     set(0,'currentfigure',figureHandle_1); hold on;
%     for k = 1:length(lines)
%         xy = [x_current(lines(k).point1(2),lines(k).point1(1)),y_current(lines(k).point1(2),lines(k).point1(1));...
%             x_current(lines(k).point2(2),lines(k).point2(1)),y_current(lines(k).point2(2),lines(k).point2(1))];
%         plot3(xy(:,1),xy(:,2),1*ones(size(xy,1),1),'LineWidth',2,'Color','green')
%     end
%
%
% end




