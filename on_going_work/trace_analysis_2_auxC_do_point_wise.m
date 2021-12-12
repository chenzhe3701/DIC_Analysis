
% do similar thing point-wise for grain 296.

iS = find(gIDwithTrace == 296); % for debugging.

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


% find vectors for cluster, using ind
ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));


% ============= clustering data.  grain-197 is a good example showing that kmeans seems to be better than gmModels =====================
% first, predict centroid

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
tLabel = (nss+1 : nss+ntwin)';         % twin system number
tSF = cPred(nss+1:nss+ntwin,2)';       % twin schmid factor
tStrain = cPred(nss+1:nss+ntwin,3:5);      % twin strain components


% ==================== continue method-1, assign cluster to twin system on the fly ==========================
twinMapLocal = zeros(size(exx_local));          % local map to record twin_system_number
sfMapLocal = zeros(size(exx_local));            % local map to record schmid_factor
disSimiMapLocal = zeros(size(exx_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
scoreMapLocal = zeros(size(exx_local));

for index = ind
    
    pdistCS = pdist2([exx_local(index), exy_local(index), eyy_local(index)], tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
    
    % score boundary y=kx+b passes [pdist2,sf] = [0, 0.15] and [0.05, 0.5]
    score = pdistCS * 7 + (0.5 - tSF);
    [m_score, ind_t] = nanmin(score,[],2);
    m_dist = pdistCS(ind_t);
    
    if isnan(m_dist)
        tsNum = [];
    else
        tsNum = tLabel(ind_t);               % match cluster to this twin system
    end
    if tsNum > nss
        % Note global assign like: twinMap(ind_c_on_map) = tsNum; is very slow.
        twinMapLocal(index) = tsNum;    % assign twinSysNum to the region in the local map
        sfMapLocal(index) = tSF(ind_t);
        disSimiMapLocal(index) = m_dist;
        scoreMapLocal(index) = m_score;
    end
end
scoreMapLocal(scoreMapLocal==0) = nan;

if 1 == debugTF
    exx_local(ID_local~=ID_current) = nan;
    exy_local(ID_local~=ID_current) = nan;
    eyy_local(ID_local~=ID_current) = nan;
    
    myplot(x_local,y_local,exx_local);
    %     set(gca,'fontsize',18);
    %     myplot(x_local,y_local,exy_local);
    %     set(gca,'fontsize',18);
    %     myplot(x_local,y_local,eyy_local);
    %     set(gca,'fontsize',18);
    myplot(x_local,y_local,clusterNumMapLocal);
    %     myplot(x_local,y_local,twinMapLocal); caxis([nss,nss+ntwin]);
    %     myplot(x_local,y_local,sfMapLocal);
    %     myplot(x_local,y_local,disSimiMapLocal);
    myplot(x_local,y_local,scoreMapLocal);
end







