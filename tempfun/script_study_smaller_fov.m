% chenzhe, 2018-01-31
% based on global EBSD-SEM control points, project a single FOV so that
% EBSD can be aligned with this SEM

%% prepare data
clear;
addChenFunction;

[fileSetting,pathSetting] = uigetfile(pwd,'select setting file which contains sampleName, stopNames, FOVs, translations, etc');
sampleName = [];    % such as 'Ti7Al_#B6'
cpEBSD = [];    % control points on EBSD image (unit um !!!!)
cpSEM = [];     % control points on SEM image (unit pixel)
sampleMaterial = [];  % such as 'Ti' or 'Mg'
stressTensor = [];
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

B = 1;   % 'B' for 'base', to handle if it's 0/1-based index.  But B=1 for 0-based. B=0 for 1-based.  When iR, iC is used with FOV, transX, ... add this B.

% data files
[EBSDfileName1, EBSDfilePath1] = uigetfile('D:\WE43_T6_C1_EBSD\*.txt','choose the EBSD file (txt format, from type-1 grain file)');
[EBSDfileName2, EBSDfilePath2] = uigetfile([EBSDfilePath1,'.txt'],'choose the EBSD file (txt format, from type-2 grain file)');
[strainFileName, strainFilePath] = uigetfile(['D:\WE43_T6_C1_insitu_compression\byFov\','*.mat'],'choose the target strain file');
[translationFileName,translationFilePath] = uigetfile(['D:\WE43_T6_C1_insitu_compression\stitched_img\translations_searched_vertical_stop_0','.mat'],'select translation');
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

[s0,s1]=regexp(strainFileName,'s\d+_');
iE = str2double(strainFileName(s0+1:s1-1));

[r0,r1]=regexp(strainFileName,'r\d+c');
[c0,c1]=regexp(strainFileName,'c\d+.');
iR = str2double(strainFileName(r0+1:r1-1));
iC = str2double(strainFileName(c0+1:c1-1));
% iR = inputdlg('input row number of fov, iR:','Input window',1,{'1'});
% iR = str2double(iR{1});
% iC = inputdlg('input col number of fov, iC:','Input window',1,{'2'});
% iC = str2double(iC{1});

%%
% Strain data
strainData = load([strainFilePath, strainFileName]);

sigma = strainData.sigma;
indnan = (sigma==-1);

try
    X = strainData.x;       % X, Y are position of SEM system.
    Y = strainData.y;
catch
    X = strainData.X;
    Y = strainData.Y;
end
try
    u = strainData.u;       % X, Y are position of SEM system.
    v = strainData.v;
catch
    u = strainData.U;
    v = strainData.V;
end
try
    exx = strainData.exx_Lagrange;
    exy = strainData.exy_Lagrange;
    eyy = strainData.eyy_Lagrange;
catch
    exx = strainData.exx;
    exy = strainData.exy;
    eyy = strainData.eyy;
end
try
    exy_corrected = strainData.exy_corrected;
catch
    exy_corrected = 0;
end
if exist('exy_corrected','var')&&(1==exy_corrected)
    disp('================= exy already corrected ! ========================');
    exy_corrected = 1;
else
    disp('================= exy being corrected here ! =======================');
    exy = -exy;
    exy_corrected = 1;
end
exx(indnan) = nan;
exy(indnan) = nan;
eyy(indnan) = nan;

% clean up, remove noise
exx(sigma==-1) = nan;
exy(sigma==-1) = nan;
eyy(sigma==-1) = nan;
qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));

remove_outlier = 1;
if remove_outlier
exx(ind_outlier) = nan;
exy(ind_outlier) = nan;
eyy(ind_outlier) = nan;
end




dic_step = Y(2)-Y(1);
% dic_step = inputdlg('input dic step size, dic_step:','Input window',1,{'1'});
% dic_step = str2double(dic_step{1});

load([translationFilePath,translationFileName]);    % load translation data
transX = dic_step * round(transX/dic_step);
transY = dic_step * round(transY/dic_step);

% Basically, this is the smallest x/y position in each FOV.
% For Vic2D-6 this should be 0.
% For Vic2D-2009, with subset=21, step = 5, you could have x=12,17,..., so
% in this case xyi=2.  Then also you should have the same x=12,17,..., for
% all your FOVs.
xyi = 0;

% This takes care of the condition when you have a translation<0.
% Since m_4 was modified to makeAllTransPositive, no need to worry here.
xi_mesh = xyi + min(transX(:));
yi_mesh = xyi + min(transY(:));

X = X + xi_mesh + transX(iR+B,iC+B);    % put local X,Y into global position, simply add diff.
Y = Y + yi_mesh + transY(iR+B,iC+B);


[EBSDdata1,EBSDheader1] = grain_file_read([EBSDfilePath1, EBSDfileName1]);
[EBSDdata2,EBSDheader2] = grain_file_read([EBSDfilePath2, EBSDfileName2]);
columnIndex1 = find_variable_column_from_grain_file_header(EBSDheader1,...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
columnIndex2 = find_variable_column_from_grain_file_header(EBSDheader2,...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

% read type-2 grain file and get average info for grains
gID = EBSDdata2(:,columnIndex2(1));
gPhi1 = EBSDdata2(:,columnIndex2(2));
gPhi = EBSDdata2(:,columnIndex2(3));
gPhi2 = EBSDdata2(:,columnIndex2(4));
gCenterX = EBSDdata2(:,columnIndex2(5));
gCenterY = EBSDdata2(:,columnIndex2(6));
gNNeighbors = EBSDdata2(:,columnIndex2(7));
gDiameter = EBSDdata2(:,columnIndex2(8));
gArea = EBSDdata2(:,columnIndex2(9));
gEdge = EBSDdata2(:,columnIndex2(10));
gNeighbors = EBSDdata2(:,(columnIndex2(7)+1):(size(EBSDdata2,2)));

% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
% Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
x = EBSDdata1(:,columnIndex1(5));
y = EBSDdata1(:,columnIndex1(6));
unique_x = unique(x(:));
ebsdStepSize = unique_x(2) - unique_x(1);
mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

phi1 = reshape(EBSDdata1(:,columnIndex1(2)),mResize,nResize)';
phi = reshape(EBSDdata1(:,columnIndex1(3)),mResize,nResize)';
phi2 = reshape(EBSDdata1(:,columnIndex1(4)),mResize,nResize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end
x = reshape(EBSDdata1(:,columnIndex1(5)),mResize,nResize)';
y = reshape(EBSDdata1(:,columnIndex1(6)),mResize,nResize)';
ID = reshape(EBSDdata1(:,columnIndex1(1)),mResize,nResize)';

clear EBSDdata1 EBSDdata2 strainData;

% align euler angle coordinates to SEM, so later on, phiSys can be set to [0 0 0]. Do this before interp so there are few points to rotate
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
% [q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));
eulerAligned = 1;

% This defines the overlay relationship, ebsdpoint(x,y) * tMatrix = sempoint(x,y)
tform = make_average_transform('projective',cpEBSD,cpSEM);
% tMatrix = tform.tdata.T;
% tInvMatrix = tform.tdata.Tinv;

% Align EBSD to SEM
% forwarded coordinates
% [x_EBSD_fwd, y_EBSD_fwd] = tformfwd(tform,x,y);

ID = interp_data(x,y,ID,X,Y,tform,'interp','nearest');
% phi1 = interp_data(x,y,phi1,X,Y,tform,'interp','nearest');
% phi = interp_data(x,y,phi,X,Y,tform,'interp','nearest');
% phi2 = interp_data(x,y,phi2,X,Y,tform,'interp','nearest');

[boundaryTF, ~, ~, ~, ~] = find_boundary_from_ID_matrix(ID);
boundaryTFB = boundaryTF;  % make it thicker
for iThick = 1:4
    boundaryTFB = grow_boundary(boundaryTFB);       % grow boundary thicker
end

save([saveDataPath,'local_data_',sampleName,'_s',num2str(iE),'_r',num2str(iR),'c',num2str(iC),'.mat']);  % record if eulerAligned.

%% Analysis
neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.

% select the target grain
ID_current = 195;   % 541

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
sigma_local = sigma(indR_min:indR_max, indC_min:indC_max);
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
x_local = X(indR_min:indR_max, indC_min:indC_max);
y_local = Y(indR_min:indR_max, indC_min:indC_max);
ID_local = ID(indR_min:indR_max, indC_min:indC_max);

% (0) can smooth
% exx_local = colfilt(exx_local, [3 3], 'sliding', @(x) nanmean(x,1));
% exy_local = colfilt(exy_local, [3 3], 'sliding', @(x) nanmean(x,1));
% eyy_local = colfilt(eyy_local, [3 3], 'sliding', @(x) nanmean(x,1));

% find vectors for cluster, using ind
ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));
exx_t = exx_local(ind);
exy_t = exy_local(ind);
eyy_t = eyy_local(ind);
data_t = [exx_t(:), exy_t(:), eyy_t(:)];
%     [data_zn,mean0,std0] = zero_normalize_column(data_t);   % seems like it's better without zero_normalize.

%% determine optimum number of clusters
rng(1); % can make consistent
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

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
    
    F = eye(2) + gamma*MN2{iss};
    epsilon = (F'*F-eye(2))/2;
    
    cPred(iss,1) = iss;                                     % ss number
    cPred(iss,2) = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
    cPred(iss,3:5) = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.  Note that 'conjugated' twin system, i.e., 19 and 22, almost always show similar components!!!
end
tLabel = (nss+1 : nss+ntwin)';         % twin system number
tSF = cPred(nss+1:nss+ntwin,2)';       % twin schmid factor
tStrain = cPred(nss+1:nss+ntwin,3:5);      % twin strain components
[~,centroid_initial_ind] = max(tSF);
centroid_initial = tStrain(centroid_initial_ind,:);
%     disp(cPred);

% ======================= kmeans, determine optimum number of clusters ====================================
maxCluster = 5;
nPoints = 8100;
ind_reduce = ~isnan(sum(data_t,2));
data_reduce = data_t(ind_reduce,:);
reduce_ratio = ceil(size(data_reduce,1)/nPoints);
data_reduce = data_reduce(1:reduce_ratio:end,:);

% compare the silhouette, by actually do kmeans on down-sampled samples.
disp(['ID=',num2str(ID_current)]);
clear wssd score_avg score_c score_min neg_score_sum;
%         score_min = -1*ones(1, maxCluster);
neg_score_sum = -inf*ones(1, maxCluster);
nRep = 3;
c0 = kmeans_pp_init(data_reduce,maxCluster,nRep,centroid_initial);
for nc = 2:maxCluster
    % nRep = 3;
    % c0 = kmeans_pp_init(data_reduce,nc,nRep,centroid_initial);
    % [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',500,'start',c0);
    [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:));   % 'correlation' distance not good.
    sil_score = silhouette(data_reduce,idx);
    %             wssd(nc) = mean(sumd);
    %             score_avg(nc) = nanmean(sil_score); % avg score for the condition of nc clusters
    for ii=1:nc
        sil_this_cluster = sil_score(idx==ii);
        %                 mean_score_cluster{nc}(ii) = mean(sil_this_cluster); % silhouette for each cluster
        neg_score_cluster{nc}(ii) = sum(sil_this_cluster(sil_this_cluster<0));
    end
    %     figure; silhouette(data_reduce,idx);
    %             score_min(nc) = min(mean_score_cluster{nc});
    neg_score_sum(nc) = sum(neg_score_cluster{nc});
end
%         [~,nCluster] = max(score_min);
[~,nCluster] = max(neg_score_sum);
disp([char(9),'nCluster=',num2str(nCluster)]);

%% (1) kmeans cluster
% nCluster = 3;       % total number of clusters

debugTF = 0;
distCI = 0.035;     % distance Criterion Initial
sfCI = 0.37;        % schmid factor Criterion Initial
shearTarget = 0.1289;
shearCI = 0.05;        % shear criterion initial
costCI = 0.05;     % sqrt(cost) criterion initial
scoreCI = 0.10;          % score criterion initial, score = 2*sqrt(cost) + abs(shearFit-shearTarget)

nRep = 3;
c0 = kmeans_pp_init(data_t,nCluster,nRep,centroid_initial);
[idx, centroid, sumd] = kmeans(data_t, nCluster, 'Distance','sqeuclidean','MaxIter',1000,'start',c0);   % 'correlation' distance not good.
clusterNumMapLocal = zeros(size(x_local));      % record raw clusterNumberMap
clusterNumMapLocal(ind) = idx;

cLabel = (1:nCluster)';     % cluster number (actually label, but labe=number)
cCen = centroid;

options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');
clear shear cost;
for iCluster=1:nCluster
    for iss = (nss+1):(nss+ntwin)
        [shear(iCluster,iss),cost(iCluster,iss)]=fminunc(@(x) sqrt(sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
            [centroid(iCluster,1),centroid(iCluster,2);centroid(iCluster,2),centroid(iCluster,3)]).^2))), 0.1298, options);
    end
end
cShear = shear(:,nss+1:nss+ntwin);       % cluster centroid's expected shear
cCost = cost(:,nss+1:nss+ntwin);      % cluster centroid fitted shear -induced cost

% =============== continue method-2, assign cluster to twin system on the fly ==============================
twinMapLocal_2 = zeros(size(exx_local));          % local map to record twin_system_number
shearMapLocal = zeros(size(exx_local));            % local map to record schmid_factor
sfMapLocal_2 = zeros(size(exx_local));       % local map to record dissimilarity between measured_strain and assigned_twin_system_theoretical_strain
costMapLocal = zeros(size(exx_local));       % local map to record cost
for iCluster = 1:nCluster
    cNum = cLabel(iCluster);
    indClusterLocal = (clusterNumMapLocal==cNum);
    
    shearFit = cShear(iCluster,:);
    shearFit(tSF < sfCI) = nan;            % [criterion-1] SF must > 0.35
    shearFit(abs(shearFit-shearTarget) > shearCI) = nan;      % [criterion-2] shear_difference should be < 0.5
    
    costFit = cCost(iCluster,:);
    costFit(costFit > costCI) = nan;     % [criterion-3] costFit should be < 0.025.  Change to 0.035 empirically
    
    score = abs(shearFit-shearTarget)+2*costFit;     % [criterion-4] the criterion is a combination of shearFit and sqrt(cost)
    score(costFit*2+abs(shearFit-shearTarget) > scoreCI) = nan;  % [criterion-5] combination of shearDiff and sqrtCost should also be limited
    
    [m_score,ind_t] = min(score,[],2);
    % [m_shear_diff,ind_t] = min(abs(shearFit-0.1289),[],2);    % if just look at the shear.
    
    if isnan(m_score)
        tsNum = [];
    else
        tsNum = tLabel(ind_t);               % match cluster to this twin system
    end
    if tsNum > nss
        twinMapLocal_2(indClusterLocal) = tsNum;
        sfMapLocal_2(indClusterLocal) = tSF(ind_t);
        shearMapLocal(indClusterLocal)  = cShear(iCluster,ind_t);
        costMapLocal(indClusterLocal) = cCost(iCluster,ind_t);
    end
end


myplot(exx_local);
myplot(clusterNumMapLocal);
myplot(twinMapLocal_2);



