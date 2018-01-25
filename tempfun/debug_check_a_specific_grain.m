
% for debug, check result for a specific grain
% select the target grain

ID_current = 668;
% ID_current = 986; 
% ID_current = 1085;
% ID_current = 1115;
% ID_current = 129;

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

% ====== this eliminates outliers , but seems not necessary for most of the grains, and actually makes them worse !!!  ==========================           
if 0
    ii = ~isnan(exx_local(:));
    ab = quantile(exx_local(ii), [0.01, 0.99]); a = ab(1); b = ab(2);
    exx_local(exx_local<a) = nan;
    exx_local(exx_local>b) = nan;
    
    ii = ~isnan(exy_local(:));
    ab = quantile(exy_local(ii), [0.01, 0.99]); a = ab(1); b = ab(2);
    exy_local(exy_local<a) = nan;
    exy_local(exy_local>b) = nan;
    
    ii = ~isnan(eyy_local(:));
    ab = quantile(eyy_local(ii), [0.01, 0.99]); a = ab(1); b = ab(2);
    eyy_local(eyy_local<a) = nan;
    eyy_local(eyy_local>b) = nan;
end
% =========================================

% find vectors for cluster, using ind
ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));
exx_t = exx_local(ind);
exy_t = exy_local(ind);
eyy_t = eyy_local(ind);
data_t = [exx_t(:), exy_t(:), eyy_t(:)];
%     [data_zn,mean0,std0] = zero_normalize_column(data_t);   % seems like it's better without zero_normalize.

myplot(x_local,y_local,exx_local);
myplot(x_local,y_local,exy_local);
myplot(x_local,y_local,eyy_local);

% ============= clustering data.  grain-197 is a good example showing that kmeans seems to be better than gmModels =====================
% (1) kmeans cluster
nCluster = 4;       % total number of clusters
[idx, centroid, sumd] = kmeans(data_t, nCluster, 'Replicates', 3, 'Distance','sqeuclidean','MaxIter',1000);   % 'correlation' distance not good.
clusterNumMapLocal = zeros(size(x_local,1),size(x_local,2));      % record raw clusterNumberMap
clusterNumMapLocal(ind) = idx;

% myplot(exx_local);myplot(exy_local);myplot(eyy_local);
figure;
colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
ind = find((ID_local==ID_current) &(~isnan(exx_local))&(~isnan(exy_local))&(~isnan(eyy_local)));
scatter3(exx_local(ind),exy_local(ind),eyy_local(ind),4/3*pi,colors(clusterNumMapLocal(ind),:));

figure;plot(quantile(exy_local(:),[0:0.01:1]),'x');

myplot(clusterNumMapLocal);

