% chenzhe, 2018-02-25
% based on codes up to date, make twin analysis into different sections.
% Here, just determine the number of clusters to use, as this is slow

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
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);

gIDwithTrace = gID(~isnan(gExx));

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
    strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
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
    % remove bad data points ----------------------------------------------------  
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
    
    % gIDwithTrace = [86, 41, 193, 194, 259, 262, 197, 296, 153, 182, 378, 451, 576]; % select some obviously twinned grain for debugging
    [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
    ss = crystal_to_cart_ss(ssa,c_a);
    %%
    clear stru struNC;
    hWaitbar = waitbar(0,'finding twin region for grains ...');
    for iS =1:length(gIDwithTrace)
        %%
        %         iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        %         iS = find(gIDwithTrace == 296); % for debugging.
        try
            stru(iS)=[];
            struNC(iS)=[];
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
        ind = find((ID_local==ID_current)); % ind = find((ID_local==ID_current)&(~isnan(exx_local))&(~isnan(exy_local))&(~isnan(eyy_local)));
        exx_t = exx_local(ind);
        exy_t = exy_local(ind);
        eyy_t = eyy_local(ind);
        data_t = [exx_t(:), exy_t(:), eyy_t(:)];
        %  [data_zn,mean0,std0] = zero_normalize_column(data_t);   % seems like it's better without zero_normalize.
        
        % export data if needed ------------------------------------------------ 
        %  save(['grain_',num2str(ID_current)],'data_t');
        
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
        
        % ======================= kmeans, determine optimum number of clusters ====================================
        maxCluster = 5;
        nPoints = 8100;
        ind_reduce = ~isnan(sum(data_t,2));
        data_reduce = data_t(ind_reduce,:);
        reduce_ratio = ceil(size(data_reduce,1)/nPoints);
        data_reduce = data_reduce(1:reduce_ratio:end,:);
        
        if(~isempty(data_reduce))
            %         % use evalclusters to evaluate the umber of clusters
            %         if isempty(data_reduce)
            %             nCluster = 4;
            %         else
            %             eva = evalclusters(data_reduce,'kmeans','silhouette','klist',2:maxCluster);
            %             nCluster = eva.OptimalK;
            %         end
            
            % compare the silhouette, by actually do kmeans on down-sampled samples.
            disp(['ID=',num2str(ID_current)]);
            clear wssd  score_avg  score_cluster_mean  score_cluster_neg_sum  score_cluster_mean_min  score_neg_sum;
            %         score_min = -1*ones(1, maxCluster);
            score_neg_sum = -inf*ones(1, maxCluster);
            nRep = 3;
            c0 = kmeans_pp_init(data_reduce,maxCluster,nRep,centroid_initial);
            for nc = 2:maxCluster
                [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:));   % 'correlation' distance not good.
                sil_score = silhouette(data_reduce,idx);
                
                % wssd(nc) = mean(sumd);
                score_avg(nc) = nanmean(sil_score); % avg score for the condition of nc clusters
                
                for ii=1:nc
                    sil_this_cluster = sil_score(idx==ii);
                    score_cluster_mean{nc}(ii) = mean(sil_this_cluster); % silhouette for each cluster
                    score_cluster_neg_sum{nc}(ii) = sum(sil_this_cluster(sil_this_cluster<0));
                end
                %             figure; silhouette(data_reduce,idx);
                score_cluster_mean_min(nc) = min(score_cluster_mean{nc});
                score_neg_sum(nc) = sum(score_cluster_neg_sum{nc});
            end
            %         [~,nCluster] = max(score_cluster_mean_min);
            [~,nCluster] = max(score_neg_sum);
            disp([char(9), 'ID= ',num2str(ID_current),', nCluster=',num2str(nCluster)]);
            
            struNC(iS).gID = ID_current;
            struNC(iS).nCluster = nCluster;
            struNC(iS).score_avg = score_avg;
            struNC(iS).score_cluster_mean = score_cluster_mean;
            struNC(iS).score_cluster_neg_sum = score_cluster_neg_sum;
            struNC(iS).score_cluster_mean_min = score_cluster_mean_min;
            struNC(iS).score_neg_sum = score_neg_sum;
            
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
            
            try
                waitbar(iS/length(gIDwithTrace), hWaitbar);
            catch
            end          
        end
    end
    
    try
        close(hWaitbar);
    catch
    end
    %%
    timeStr = datestr(now,'yyyymmdd_HHMM')
    nameNC = [sampleName,'_s',num2str(STOP{iE+B}),'_nClusters_',timeStr,'.mat'];
    save([saveDataPath,nameNC],'struNC');
    disp(['finished strain: ',num2str(STOP{iE+B})]);
    
end




