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

%% Can load strain data for a specific strain level
iE = 5;
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
load(strainFile,'exx','exy','eyy');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------

%%
% for debug, create a few maps to record the criterion.
twinMap = zeros(size(exx));
sfMap = zeros(size(exx));
disSimiMap = zeros(size(exx));
twinMap_2 = zeros(size(exx));
shearMap = zeros(size(exx));
sfMap_2 = zeros(size(exx));

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

% gIDwithTrace = [86, 41, 193, 194, 259, 262, 197, 296, 153, 182, 378, 451, 576]; % select some obviously twinned grain for debugging
g_aoi_range = [];

hWaitbar = waitbar(0,'finding twin region for grains ...');
for iS =1:length(gIDwithTrace)
    close all;
    % select the target grain
    ID_current = gIDwithTrace(iS);  % id=262 for an example for WE43-T6-C1
    %     ID_current = 668;
    
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
    for iss=19:24   % for Mg
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
    %     disp(cPred);
    
    % =================== find measured strain: =======================
    ind_current = find(ID_current == gID);    % an index of row
    phi1_local = gPhi1(ind_current);
    phi_local = gPhi(ind_current);
    phi2_local = gPhi2(ind_current);
    
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
    
    %     myplot(x_local,y_local,exx_local);
    %     myplot(x_local,y_local,exy_local);
    %     myplot(x_local,y_local,eyy_local);
    
    % ============= clustering data.  grain-197 is a good example showing that kmeans seems to be better than gmModels =====================
    % (1) kmeans cluster
    nCluster = 4;       % total number of clusters
    [idx, centroid, sumd] = kmeans(data_t, nCluster, 'Distance','sqeuclidean','MaxIter',1000);   % 'correlation' distance not good.
    clusterNumMapRaw = zeros(size(x_local,1),size(x_local,2));      % record raw clusterNumberMap
    clusterNumMapRaw(ind) = idx;
    
    % ================ method-1, from theoretical twin shear --> to predict strain components as cluster center ===================
    twinMapLocal = zeros(size(x_local,1),size(x_local,2));          % local map to record twin_system_number
    sfMapLocal = zeros(size(x_local,1),size(x_local,2));            % local map to record schmid_factor
    disSimiMapLocal = zeros(size(x_local,1),size(x_local,2));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
    %     fh1 = myplot(x_local,y_local,cNMap);
    pdistCS = pdist2(centroid, cPred(:,3:5));       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
    %     pdistCS(pdistCS > 0.025) = nan;                 % [criterion-1] can do this: if a cluster center-slip system center distance is too large, this cluster shouldn't be a twin system
    [m_dist, matchingTS] = nanmin(pdistCS,[],2);    % [criterion-2] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
    matchingTS(isnan(m_dist)) = 0;
    
    % ==== change cluster number into twin system number, or 0
    for iCluster = 1:nCluster
        tsNum = matchingTS(iCluster);               % twin system number currently considered 
        if tsNum > 18
            indTwinLocal = (clusterNumMapRaw==iCluster);
            twinMapLocal(indTwinLocal) = tsNum;    % assign twinSysNum to the region in the local map
            sfMapLocal(indTwinLocal) = cPred(tsNum,2);
            disSimiMapLocal(indTwinLocal) = m_dist(iCluster);
        end
    end
    %     fh2 = myplot(x_local,y_local,cNMap);
    
    % copy identified twin system number to twinMap
    twinMap(indR_min:indR_max, indC_min:indC_max) = twinMap(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
    sfMap(indR_min:indR_max, indC_min:indC_max) = sfMap(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
    disSimiMap(indR_min:indR_max, indC_min:indC_max) = disSimiMap(indR_min:indR_max, indC_min:indC_max) + disSimiMapLocal;
    
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
    
    
    % ============ method-2, from cluster_centroid to estimated_shear, and compare if this is similar enough to theoretical shear of 0.1289 
    twinMapLocal = zeros(size(x_local,1),size(x_local,2));          % local map to record twin_system_number
    shearMapLocal = zeros(size(x_local,1),size(x_local,2));            % local map to record schmid_factor
    sfMapLocal = zeros(size(x_local,1),size(x_local,2));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
    
    options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');

    matchingTS = zeros(nCluster,1);
    s = ones(nCluster,1)*100;   % initialize with very large value.
    cost = zeros(nCluster,1);
    for iCluster=1:nCluster
        % fit the shear 's' for each cluster, by find the most likely twin system: 
        %   for each cluster
        %       for each if possible ts
        %           calculate fitted twin shear
        %           if abs(fitted twin shear - 0.1298) is smaller
        %               update cluster twin type, twin shear, cost
        for iss = 19:24
           if ~isnan(cPred(iss,1))
              [s_i,cost_i]=fminunc(@(x) sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
                  [centroid(iCluster,1),centroid(iCluster,2);centroid(iCluster,2),centroid(iCluster,3)]).^2)), 0.1298, options);
              if abs(s_i-0.1298) < abs(s(iCluster)-0.1298)
                  matchingTS(iCluster) = iss;
                  s(iCluster) = s_i;
                  cost(iCluster) = cost_i;
              end              
           end
        end
        % ==== similarly, change cluster number into twin system number, or 0
        tsNum = matchingTS(iCluster);               % twin system number currently considered
        if tsNum > -1   % this should be 18, but set to -1 for debug
            indTwinLocal = (clusterNumMapRaw==iCluster);
            twinMapLocal(indTwinLocal) = tsNum;    % assign twinSysNum to the region in the local map
            shearMapLocal(indTwinLocal) = s(iCluster);
            if tsNum>0
                sfMapLocal(indTwinLocal) = cPred(tsNum,2);
            end
        end
    end
    % copy identified twin system number to twinMap
    twinMap_2(indR_min:indR_max, indC_min:indC_max) = twinMap_2(indR_min:indR_max, indC_min:indC_max) + twinMapLocal;
    shearMap(indR_min:indR_max, indC_min:indC_max) = shearMap(indR_min:indR_max, indC_min:indC_max) + shearMapLocal;
    sfMap_2(indR_min:indR_max, indC_min:indC_max) = sfMap_2(indR_min:indR_max, indC_min:indC_max) + sfMapLocal;
    
    %     input('press to continue');
    waitbar(iS/length(gIDwithTrace), hWaitbar);
end
close(hWaitbar);
disp('start saving tempWS');
save('tempWS.mat');
disp('finished saving tempWS');

%% temp code for plot and investigate the results
myplot(exx,boundaryTFB);
myplot(disSimiMap,boundaryTFB);
myplot(twinMap,boundaryTFB); caxis([19,24]);
myplot(sfMap,boundaryTFB);

myplot(shearMap, boundaryTFB);
myplot(twinMap_2, boundaryTFB); caxis([19,24]);
myplot(sfMap_2, boundaryTFB);


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






















