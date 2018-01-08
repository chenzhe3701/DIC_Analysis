%% This part is definition
% EBSD scan step size (in um)
clc;
addChenFunction;
step_EBSD = 1;
step_dream = 2;
AVERAGED = 1;               % use averaged strain
averaged_filter_size = 601; % and the filter size to average the strain

% for reference, these are some of the 
sampleName = '';
cpEBSD = []; % control points on EBSD image (unit um)
cpSEM = [];  % control points on SEM image (unit pixel)
DIC_resolution = [];
stressTensor = [];    % stress state in Global System, be careful of its direction.  Here tensile axis is horizontal.
material = [];  % might be 'Ti'
cpStruct_dream_dic = [];    % control points which controls alignment.
cpStruct_ebsd_dic = [];
cpStruct_dream_ebsd = [];
[f,p] = uigetfile('d:\p\m\Ti7Al_#B6_setting.mat','select settings');
load_settings([p,f],'sampleName','cpEBSD','cpSEM','DIC_resolution','stressTensor','material','cpStruct_dream_dic','cpStruct_ebsd_dic','cpStruct_dream_ebsd');

cpSEM_um = cpSEM*DIC_resolution;
tform1 = maketform('projective',cpEBSD,cpSEM);   % SEM in pixel
tform2 = maketform('projective',cpEBSD,cpSEM_um);   % SEM in micron

% data files
[EBSDfileName1, EBSDfilePath1] =...
    uigetfile('D:\UMich Folder\11 ONR project\31 Experiemt data\[Analysis] Ti7Al #B6 EBSD\Ti7Al#B6_Left_GrainFile_Type_1.csv', ...
    'choose the EBSD file (csv format, from type-1 grain file)');
[EBSDfileName2, EBSDfilePath2] =...
    uigetfile('D:\UMich Folder\11 ONR project\31 Experiemt data\[Analysis] Ti7Al #B6 EBSD\Ti7Al#B6_Left_GrainFile_Type_2.csv', ...
    'choose the EBSD file (csv format, from type-2 grain file)');

[strainFileName, strainFilePath] =...
    uigetfile('F:\Ti7Al #B6 DIC_FOVs_merged_highRes_WithNaN\Ti7Al_#B6_stop_015.mat',...
    'choose the strain file (mat format)');
[avgedStrainFileName, avgedStrainFilePath] =...
    uigetfile('F:\Ti7Al #B6 DIC_FOVs_merged_highRes_Avg_VaryingFilterSize\Ti7Al_#B6_stop_015.mat',...
    'choose the strain file (mat format)');
dicFolder =...
    uigetdir('F:\Ti7Al #B6 DIC_FOVs_merged_highRes_Avg_VaryingFilterSize',...
    'choose the whole DIC folder');
dicFiles = dir([dicFolder,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

[dreamFileName, dreamFilePath] =...
    uigetfile('D:\3D data\MyOutputData\Ti7Al_B6_104_Layer_step_3a.dream3d',...
    'choose dream3D file');

%  define slip systems, format: [slip PLANE, slip DIRECTION], and change to unit vector
[ssa, c_a] = define_SS(material,'twin');   % slip system for Ti
nSS = size(ssa,3);
% for Euler angle related analysis:
phiSys = [0, 0, 0];  % angle to rotate sample system to align with euler/EBSD system.  Because I've decided to align euler with sample during data importing, this is now set to [0 0 0] !!!
phiError = [0, 0, 0];  % angle error caused by difference in placing the sample when doing experiment & when doing EBSD.
% i.e., We need to rotate EBSD map into In-situ SEM image orientation.
% This phiError is the euler angle to rotate SEM position into the EBSD position.
% Note: If you only need to rotate SEM image CCW to match EBSD map, then the angle phiError(1) is Positive.

save([sampleName,'_WS1.mat']);

%% Load DIC data
loaded = load([strainFilePath,strainFileName]);
% sigma = loaded.sigma;
% u = loaded.u;
% v = loaded.v;
x = loaded.x;
y = loaded.y;

if AVERAGED == 1   % directly load averaged
    loaded = load([avgedStrainFilePath,avgedStrainFileName]);
    exx = loaded.exxA;
    exy = loaded.exyA;
    eyy = loaded.eyyA;
else
    exx = loaded.exx;
    exy = loaded.exy;
    eyy = loaded.eyy;
end
clear loaded;

%%
% load files, EBSD type-2
EBSDdata2 = csvread([EBSDfilePath2, EBSDfileName2],1,0);
columnIndex2 = find_variable_column_from_CSV_grain_file(...
    EBSDfilePath2, EBSDfileName2,...
    {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

for iSS = 1:1:nSS            % Change n & m to unit vector
    planeNormal = [ssa(1,1,iSS)  (ssa(1,2,iSS)*2+ssa(1,1,iSS))/sqrt(3)  ssa(1,4,iSS)/c_a]; % Plane normal /c_a, into a Cartesian coordinate
    slipDirection = [ssa(2,1,iSS)*3/2  (ssa(2,1,iSS)+ssa(2,2,iSS)*2)*sqrt(3)/2  ssa(2,4,iSS)*c_a]; % Slip direction *c_a, into a Cartesian coordinate
    ss(1,:,iSS) = planeNormal/norm(planeNormal);  % slip plane, normalized
    ss(2,:,iSS) = slipDirection/norm(slipDirection);  % slip direction, normalized
end
clear iSS;

gID = EBSDdata2(:,columnIndex2(1));
gPhi1 = EBSDdata2(:,columnIndex2(2));
gPhi = EBSDdata2(:,columnIndex2(3));
gPhi2 = EBSDdata2(:,columnIndex2(4));
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', 90,180,0); % align euler angle to sample reference frame -------------- align. Umich data is actually setting-1 !! 

gCenterX = EBSDdata2(:,columnIndex2(5));    % this can be modified due to projection
gCenterY = EBSDdata2(:,columnIndex2(6));    % this can be modified due to projection
gNNeighbors = EBSDdata2(:,columnIndex2(7));
gDiameter = EBSDdata2(:,columnIndex2(8));   % this can be modified due to projection
gArea = EBSDdata2(:,columnIndex2(9));       % this can be recalculated due to projection
gEdge = EBSDdata2(:,columnIndex2(10));
gNeighbors = EBSDdata2(:,(columnIndex2(7)+1):(size(EBSDdata2,2)));

% construct grain neighbor structure
neighborStruct = construct_neighbor_structure(...
    EBSDfilePath2,EBSDfileName2);
% misorientationStruct.g1 = [1;2;3], misorientationStruct.g2{i} = [2;3;...],
% misorientationStruct.misorientation{i}=[5d;75d;...]
misorientationStruct = construct_misorientation_structure(...
    neighborStruct, gPhi1, gPhi, gPhi2);

% find average schmidFactor for grain
[gSFBa,gSFPr,gSFPy,gSFPyCA,gSFTT,...
    gSFBaNum,gSFPrNum,gSFPyNum,gSFPyCANum,gSFTTNum]...
    = calculate_SFs(gPhi1, gPhi, gPhi2,...
    ss, phiSys(1), phiSys(2), phiSys(3), phiError(1), phiError(2), phiError(3), stressTensor, material);
clear gSFBaNum gSFPrNum gSFPyNum gSFPyCANum gSFTTNum;
%% EBSD data, from type-1 grain file. (column, data) pair:
% (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
EBSDdata1 = csvread([EBSDfilePath1, EBSDfileName1],1,0);
columnIndex1 = find_variable_column_from_CSV_grain_file(...
    EBSDfilePath1, EBSDfileName1,...
    {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
% read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
X = EBSDdata1(:,columnIndex1(5));
Y = EBSDdata1(:,columnIndex1(6));
mResize = (max(X(:)) - min(X(:)))/step_EBSD + 1;
nResize = (max(Y(:)) - min(Y(:)))/step_EBSD + 1;

phi1 = reshape(EBSDdata1(:,columnIndex1(2)),mResize,nResize)';
phi = reshape(EBSDdata1(:,columnIndex1(3)),mResize,nResize)';
phi2 = reshape(EBSDdata1(:,columnIndex1(4)),mResize,nResize)';
% change it to degrees, if necessary
if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
    phi1 = phi1*180/pi();
    phi = phi*180/pi();
    phi2 = phi2* 180/pi();
end

[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', 90,180,0); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
[q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles

X = reshape(EBSDdata1(:,columnIndex1(5)),mResize,nResize)';
Y = reshape(EBSDdata1(:,columnIndex1(6)),mResize,nResize)';
ID = reshape(EBSDdata1(:,columnIndex1(1)),mResize,nResize)';
edge = reshape(EBSDdata1(:,columnIndex1(7)),mResize,nResize)';

clear columnIndex1 mResize nResize;
%% align / select data
ID_for_DIC = interp_data(X,Y,ID,x,y,tform1,'interp','nearest');
xum_for_DIC = x * DIC_resolution;
yum_for_DIC = y * DIC_resolution;
if AVERAGED ~= 1
    [exx, exy, eyy] = average_by_varying_filter(exx, exy, eyy, ID_for_DIC, averaged_filter_size);
end

% select a portion of the data, lower resolution/amount of DIC data,
% because the EBSD only has limited resolution/amount of data.
step_DIC_target = 1;    % in micron

[xum_DIC, yum_DIC] = meshgrid(xum_for_DIC(1,1):step_DIC_target:xum_for_DIC(1,end),...
    yum_for_DIC(1,1):step_DIC_target:yum_for_DIC(end,1));

ID_DIC = interp_data(xum_for_DIC, yum_for_DIC, ID_for_DIC, xum_DIC, yum_DIC, [], 'interp', 'nearest');

filterSize = round(step_DIC_target/DIC_resolution);
filter = fspecial('average',filterSize);
exx_DIC = filter2(filter,exx,'same');
exy_DIC = filter2(filter,exy,'same');
eyy_DIC = filter2(filter,eyy,'same');

exx_DIC = interp_data(xum_for_DIC,yum_for_DIC,exx_DIC,xum_DIC,yum_DIC,[],'interp','nearest');
exy_DIC = interp_data(xum_for_DIC,yum_for_DIC,exy_DIC,xum_DIC,yum_DIC,[],'interp','nearest');
eyy_DIC = interp_data(xum_for_DIC,yum_for_DIC,eyy_DIC,xum_DIC,yum_DIC,[],'interp','nearest');
e_DIC =  sqrt(2/3*(exx_DIC.^2+eyy_DIC.^2+2*exy_DIC.^2));    % Equivalent Von Mises Strain
% note:
% e3 = sqrt(0.25*(exx-eyy).^2+exy.^2);
% e1 = (exx+eyy)/2 + e3;
% e2 = (exx+eyy)/2 - e3;
clear EBSDdata1 EBSDdata2 ID_for_DIC xum_for_DIC yum_for_DIC filterSize filter;
save([sampleName,'_WS2.mat']);
%% for pre-EBSD data
% find if a point is on grain boundary (boundaryTF).
% If it is a grain boundary, then which grain it belongs to (boundaryID),
% which grain it is adjacent to (neighborID)
[boundaryTF, boundaryID, neighborIDonGb, tripleTF, tripleID]...
    = find_boundary_from_ID_matrix(ID);

% For each point, calculate its distance to its grain boundary (unit is the same as x,y,gCenterX,gCenterY)
% In addition, find its nearest neighboring grain (i.e., make neighborID complete)
[dToBoundary,dToCenter,dToTriple, rToBoundary,rToCenter,rToTriple, neighborID]...
    = find_d_to_boundary_center_triple(X,Y,ID,...
        boundaryID,tripleID,neighborIDonGb,gID,...
        gCenterX,gCenterY,gDiameter);

[sfBa,sfPr,sfPy,sfPyCA,sfTT,...
    sfBaNum,sfPrNum,sfPyNum,sfPyCANum,sfTTNum]...
    = calculate_SFs(phi1, phi, phi2, ss, phiSys(1),phiSys(2),phiSys(3), phiError(1),phiError(2),phiError(3), stressTensor, material);
clear sfBaNum sfPrNum sfPyNum sfPyCANum sfTTNum;    % no need

% for DIC data
% find if a point is on grain boundary (boundaryTF).
% If it is a grain boundary, then which grain it belongs to (boundaryID),
% which grain it is adjacent to (neighborID)
[boundaryTF_DIC, boundaryID_DIC, neighborIDonGb_DIC, tripleTF_DIC, tripleID_DIC]...
    = find_boundary_from_ID_matrix(ID_DIC);

%% Dream3D data
% cell/pixel/voxel data
ID_dream = load_dream_data([dreamFilePath, dreamFileName],'FeatureIds','CellData',1);

[boundaryTF_dream, boundaryID_dream, neighborIDonGb_dream, tripleTF_dream, tripleID_dream]...
    = find_boundary_from_ID_matrix(ID_dream);

[phi1_dream,phi_dream,phi2_dream] = load_dream_data([dreamFilePath, dreamFileName],'EulerAngles','CellData',1);
phi1_dream = phi1_dream*180/pi();
phi_dream = phi_dream*180/pi();
phi2_dream = phi2_dream* 180/pi();
[q0_dream,q1_dream,q2_dream,q3_dream,phi1_dream,phi_dream,phi2_dream] = regulate_euler_quat(phi1_dream,phi_dream,phi2_dream);

[gPhi1_dream,gPhi_dream,gPhi2_dream] = load_dream_data([dreamFilePath, dreamFileName],'gAvgEuler','CellFeatureData');
gPhi1_dream = gPhi1_dream*180/pi();
gPhi_dream = gPhi_dream*180/pi();
gPhi2_dream = gPhi2_dream* 180/pi();

gID_dream = (0:length(gPhi1_dream)-1)';     % the gID in dream3D file has redundant grain numbers. !!!

nCells_dream = load_dream_data([dreamFilePath,dreamFileName],'numCells','CellData',1);
gNCells_dream = load_dream_data([dreamFilePath,dreamFileName],'gNumCells','CellFeatureData');

% temporarily treat like this, using step=2. Cannot figure out the meaning of 'X position', 'Y position' in dream3d data. 
[X_dream, Y_dream] = meshgrid(0:step_dream:(size(ID_dream,2)-1)*2, 0:step_dream:(size(ID_dream,1)-1)*2);

%% align data
s_ebsd = struct('ID',ID,'gID',gID,'gEuler',[gPhi1,gPhi,gPhi2]);
s_dic = struct('ID',ID_DIC,'gID',gID,'gEuler',[gPhi1,gPhi,gPhi2]);
s_dream = struct('ID',ID_dream,'gID',(0:size(gPhi1_dream,1)-1)','gEuler',[gPhi1_dream,gPhi_dream,gPhi2_dream]);     % dream3d's ID starts from 0 !!!


%%
% align EBSD and DIC.
[indMatch_dream_dic, cpStruct_dream_dic] = align_IDs(s_dream, s_dic, cpStruct_dream_dic);
[indMatch_ebsd_dic, cpStruct_ebsd_dic] = align_IDs(s_ebsd, s_dic, cpStruct_ebsd_dic);
[indMatch_dream_ebsd, cpStruct_dream_ebsd] = align_IDs(s_dream, s_ebsd, cpStruct_dream_ebsd);
indexMatch = combine_match(indMatch_ebsd_dic,fliplr(indMatch_dream_dic));    
    
save([sampleName,'_WS3.mat']);
%%  (1) look at grain diameter 2d vs 3d

diameter = assign_field_to_cell(ID_DIC,gID,gDiameter);
diameter_3D = 2*(3/4/pi*nCells_dream*step_dream^3).^(1/3);

diameterB = zeros(size(diameter));
for ii=1:length(diameterB(:))
    ind_M = indexMatch(ii,2);
    ind_B = indexMatch(ii,3);
    if (ind_M>0)&&(ind_B>0)
        diameterB(ind_M) = diameter_3D(ind_B);
    end
end

%%  (2) try load dream3d multiple layers, compare local misorientation

nCells_dream = load_dream_data([dreamFilePath,dreamFileName],'numCells','CellData',1);

%% calculate local orientation spread
los = calculate_los(q0,q1,q2,q3,20,30);








