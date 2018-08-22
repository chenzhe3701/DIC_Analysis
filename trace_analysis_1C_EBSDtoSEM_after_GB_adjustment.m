
% chenzhe, 2018-07-28
% After adjusting the grain boundary overlap, it is better to re-align EBSD
% to SEM

clear;
clc;
addChenFunction;

[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
sampleName = [];    % such as 'Ti7Al_#B6'
cpEBSD = [];    % control points on EBSD image (unit um !!!!)
cpSEM = [];     % control points on SEM image (unit pixel)
sampleMaterial = [];  % such as 'Ti' or 'Mg'
stressTensor = [];
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% data files
[EBSDfileName1, EBSDfilePath1] = uigetfile('*.txt','choose the EBSD file (txt format, from type-1 grain file)');
[EBSDfileName2, EBSDfilePath2] = uigetfile([EBSDfilePath1,'.txt'],'choose the EBSD file (txt format, from type-2 grain file)');
[strainFileName, strainFilePath] = uigetfile([EBSDfilePath1,'.mat'],'choose one of the strain file (mat format) for aligning purpose');
[gbReAlignedFileName, gbReAlignedFilePath] = uigetfile(['D:\p\m\DIC_Analysis\setting_for_real_samples\',sampleName,'_boundary_model.mat'],'choose the file of the re-aligned grain boundary data');

saveDataPath = [uigetdir(pathSetting,'choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

umPerPtDIC = 360/4096*5;
min_valid_gID = 2;  % e.g., grain # 1 represents the frame of the EBSD area, so it is not a valid grain.

try
    save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],'-append');
catch
    save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
end

demo = 1;   % some transformations can be demonstrated as not the best choice. Use demo=1 to plot the results to confirm.
%% Read raw EBSD data and DIC data.  The alignment will be performed in next section.
% Read EBSD grain file. '.txt' format for now.
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
ID_0 = ID;  % keep a copy, maybe useful
edge = reshape(EBSDdata1(:,columnIndex1(7)),mResize,nResize)';


% Strain data
strainData = load([strainFilePath, strainFileName]);
try
    X = strainData.x;       % X, Y are position of SEM system.
    Y = strainData.y;
catch
    X = strainData.X;
    Y = strainData.Y;
end
sigma = strainData.sigma;
try
    exx = strainData.exx_Lagrange;
    exy = strainData.exy_Lagrange;
    eyy = strainData.eyy_Lagrange;
catch
    exx = strainData.exx;
    exy = strainData.exy;
    eyy = strainData.eyy;
end

% When load new strain data, First clear 'exy_corrected', then load
% 'exy_corrected'.  If there is no such variable, it leads to a warning
% rather than error.  Then you know if it was corrected, and can generate
% and save the correct 'exy_corrected' variable.
clear exy_corrected;
load([strainFilePath, strainFileName], 'exy_corrected');

if exist('exy_corrected','var')&&(1==exy_corrected)
    disp('================= exy already corrected ! ========================');
    exy_corrected = 1;
else
    disp('================= exy being corrected here ! =======================');
    exy = -exy;
    exy_corrected = 1;
end
clear EBSDdata1 EBSDdata2 strainData;

load([sampleName,'_boundary_model.mat'], 'ID_aligned', 'ID_link_additional');   % this is the ID that was adjusted to align well with the strain map.


%% Use the selected adjusted grains to define an averaged transformation for alignment, ebsdpoint(x,y) * tMatrix = sempoint(x,y)
ids_1 = unique(ID_aligned); % all ids
ids_2 = [unique(ID_aligned(1,:)), unique(ID_aligned(end,:))];
ids_3 = [unique(ID_aligned(:,1)); unique(ID_aligned(:,end))];
ids_4 = unique([ids_2(:);ids_3(:)]);    % ids on the border
ids_5 = unique(ID_link_additional(:));    % ids that are related to the new/splitted old ids

ids_for_avg_tform = ids_1(~ismember(ids_1,[ids_4(:);ids_5(:)]));    % ids used to re-align
ids_interior = ids_1(~ismember(ids_1,ids_4));    % adjusted ids that are not on the edge

% (0) extract some statistics for the grains in the adjusted AOI
[gID_AOI,~,IC] = unique(ID_aligned); % here, unique_ID_AOI contains unique IDs on the gb_adjusted map in the AOI
accumNumPts = accumarray(IC,1);
accumX = accumarray(IC,X(:));
accumY = accumarray(IC,Y(:));

gCenterX_AOI = accumX./accumNumPts;     % recalculate grain center directly from map
gCenterY_AOI = accumY./accumNumPts;

gArea_AOI = accumNumPts * umPerPtDIC * umPerPtDIC;
gDiameter_AOI = (4*gArea_AOI/pi).^(0.5);

% (1) find out the grain center on the DIC map
ind = ismember(gID_AOI,ids_for_avg_tform); % because both 'unique_ID_AOI' and 'ids_for_avg_tform' are sorted
cpSEM = [gCenterX_AOI(ind),gCenterY_AOI(ind)];

% (2) find out the grain center on the old EBSD map
ind = ismember(gID,ids_for_avg_tform);    % because both 'gID' and 'ids_for_avg_tform' are sorted
cpEBSD = [gCenterX(ind),gCenterY(ind)];

% (3) make average tform
tform = make_average_transform('affine',cpEBSD,cpSEM);
tMatrix = tform.tdata.T;
tInvMatrix = tform.tdata.Tinv;

%% align euler angle coordinates to SEM, so later on, phiSys can be set to [0 0 0]. Do this before interp so there are few points to rotate
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
% [q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));

eulerAligned = 1;
save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],'eulerAligned','-append');  % record if eulerAligned.

%% Align EBSD to SEM
% ID_aligned is directly available
boundaryTF_aligned = find_boundary_from_ID_matrix(ID_aligned);

% Roughly, can used an averaged transformation again, based on the adjusted grain boundary map
% For grains that are in the AOI, we could directly found from ID_realigned map later.
[gCenterX, gCenterY] = tformfwd(tform,[gCenterX, gCenterY]);

% We only need to re-align Euler angles, to put on top of re-aligned data
% However, it can be shown that after average transform again, usually the average transformed data still do not align with the adjusted boundary well.
% So, point-wise data need to be better transformed, possibly grain-by-grain. Otherwise, maybe just use the average data.
% Therefore, this part is actually not needed.
if demo==1
    phi1_aligned = interp_data(x,y,phi1,X,Y,tform,'interp','nearest');
    phi_aligned = interp_data(x,y,phi,X,Y,tform,'interp','nearest');
    phi2_aligned = interp_data(x,y,phi2,X,Y,tform,'interp','nearest');
    ID_avg_tformed = interp_data(x,y,ID_0,X,Y,tform,'interp','nearest');
    myplot(rem(ID_avg_tformed,5),boundaryTF_aligned);
end

%% after this, things in the EBSD data that need to be modified: ID, gID, gPhi1, gPhi, gPhi2, gCenterX, gCenterY, gNNeighbors, gDiameter, gArea, gEdge, gNeighbors
% make a new gID_aligned, which is the old gID + newly assigned
gID_additional = ID_link_additional(:,1);
for ii = 1:size(ID_link_additional,1)
    ind = gID==ID_link_additional(ii,2);
    gPhi1_additional(ii,1) = gPhi1(ind);
    gPhi_additional(ii,1) = gPhi(ind);
    gPhi2_additional(ii,1) = gPhi2(ind);
end

% (1) for gID and gEuler, append additional to the end
gID_aligned = [gID; gID_additional];
gPhi1_aligned = [gPhi1; gPhi1_additional];
gPhi_aligned = [gPhi;gPhi_additional];
gPhi2_aligned = [gPhi2;gPhi2_additional];

% (2) for gEdge, just directly analyze the position. If outside of AOI, label as on edge
ind = ismember(gID_aligned, ids_interior);      % index indicating interior grains
gEdge_aligned = ones(size(gID_aligned));
gEdge_aligned(ind) = 0;

% (3) for gCenterX/Y, gArea, gDiameter. In previous section, scaled old ones. Now modify the new ones in AOI
for ii = 1:size(gID_aligned,1)
    ind = find(gID_aligned(ii) == gID_AOI);
    % if in DIC area, find grain center from new map. Else, use interpolated.
    if ~isempty(ind)
        gCenterX_aligned(ii,1) = gCenterX_AOI(ind);
        gCenterY_aligned(ii,1) = gCenterY_AOI(ind);
        gArea_aligned(ii,1) = gArea_AOI(ind);
        gDiameter_aligned(ii,1) = gDiameter_AOI(ind);
    else
        ind = find(gID_aligned(ii) == gID);
        gCenterX_aligned(ii,1) = gCenterX(ii);
        gCenterY_aligned(ii,1) = gCenterY(ii);
        gArea_aligned(ii,1) = gArea(ind);
        gDiameter_aligned(ii,1) = gDiameter(ind);
    end
    
end

% show how the new gCenter differ from old ones
if demo == 1
    figure; hold on;
    plot(gCenterX,gCenterY,'.r');
    plot(gCenterX_aligned,gCenterY_aligned,'ob');
    plot(X(boundaryTF_aligned(:)==1),Y(boundaryTF_aligned(:)==1),'.k');
    legend('Interpolated','ReCalculated');
    set(gca,'ydir','reverse')
end

% (4) for gNNeighbors, gNeighbors, neighborStruct, misorientationStruct, directly find on adjusted/re-aglinged map   
[boundaryTF_aligned, boundaryID_aligned, neighborID_aligned, ~, ~,] = find_one_boundary_from_ID_matrix(ID_aligned);
pairs = unique([ID_aligned(boundaryTF_aligned==1), neighborID_aligned(boundaryTF_aligned==1)], 'rows');
ind = pairs(:,1)<min_valid_gID | pairs(:,2)<min_valid_gID;
pairs(ind,:) = [];
pairs = unique([pairs; fliplr(pairs)], 'rows');

clear gNNeighbors_aligned gNeighbors_aligned neighborStruct_aligned;

% [neighborStruct] construct grain neighbor structure S.g1 = [1;2;3], S.g2{ind_in_g1} = [2;3;...]
% [misorientationStruct] Note this is currently only for HCP!
% misorientationStruct.g1 = [1;2;3], misorientationStruct.g2{i} = [2;3;...],
% misorientationStruct.misorientation{i}=[5d;75d;...]
for ii = 1:length(gID_aligned)
    neighborStruct_aligned.g1(ii) = gID_aligned(ii);
    
    ind = gID_aligned(ii) == pairs(:,1);
    if ~isempty(ind)
        gNNeighbors_aligned(ii,1) = sum(ind);
    else
        gNNeighbors_aligned(ii,1) = 0;
    end
end

gNeighbors_aligned = zeros(length(gID_aligned), max(gNNeighbors));
for ii = 1:length(gID_aligned)
    ind = gID_aligned(ii) == pairs(:,1);
    if sum(ind)>0
        gNeighbors_aligned(ii, 1:sum(ind)) = reshape(pairs(ind,2),1,[]);
        neighborStruct_aligned.g2{ii} = reshape(pairs(ind,2),1,[]);
    else
        neighborStruct_aligned.g2{ii} = 0;
    end
end

misorientationStruct_aligned = construct_misorientation_structure(neighborStruct_aligned, gPhi1_aligned, gPhi_aligned, gPhi2_aligned);


%% Save the data
[gQ0,gQ1,gQ2,gQ3,gPhi1,gPhi,gPhi2] = regulate_euler_quat(gPhi1_aligned,gPhi_aligned,gPhi2_aligned);   % regulate the angles
% Grain average strains should be done after SEM and EBSD is aligned.
% If > 20% data points are valid in this grain, then there will be an avg value.  Otherwise, the value is NaN. ----------------- can modify.
[gExx_aligned,~] = generate_grain_avg_data(ID_aligned,gID_aligned,exx, 0.2, sigma);

disp('saving ...');
ID = ID_aligned; 
boundaryTF = boundaryTF_aligned;
boundaryTFB = boundaryTF;
for iThick = 1:4
    boundaryTFB = grow_boundary(boundaryTFB);       % grow boundary thicker
end
neighborStruct = neighborStruct_aligned;
misorientationStruct = misorientationStruct_aligned;

gArea = gArea_aligned;
gCenterX = gCenterX_aligned;
gCenterY = gCenterY_aligned;
gDiameter = gDiameter_aligned;
gExx = gExx_aligned;
gID = gID_aligned;
gNNeighbors = gNNeighbors_aligned;
gNeighbors = gNeighbors_aligned;


save([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'], 'ID','X','Y','x','y','boundaryTF','boundaryTFB','eulerAligned',...
    'gPhi1','gPhi','gPhi2','gQ0','gQ1','gQ2','gQ3',...
    'neighborStruct','misorientationStruct',...    
    'gArea','gCenterX','gCenterY','gDiameter','gExx','gID','gNNeighbors',...
    'gNeighbors','gPhi1','gPhi','gPhi2','exx','exy','eyy','sigma',...
    'ebsdStepSize','fileSetting','pathSetting',...
    'sampleName','sampleMaterial','stressTensor');




