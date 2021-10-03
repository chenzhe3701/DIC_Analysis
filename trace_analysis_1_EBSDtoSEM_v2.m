% ZheChen, 2015-11-29
%
% To performe trace analysis, the trace direction needs to be as accurate
% as possible. When performing EBSD, there might be a few degrees error in
% angle of crystal oriengtation, but there might be +10% in height scanned.
%
% Take an extreme condition, sample tilt is 15 degree off (~90 instead of
% 75), the measured orientation is 15 degrees off.  The predicted trace
% direction will also be off, but probably will be < 15 degrees.
%
% On the other hand, all traces observed on SEM need to shrink in height
% direction to align with EBSD, this means all observed traces, after
% projected on EBSD, will look like almost horizontal, that's ~90 degree
% error maximum
%
% Therefore, when performing trace analysis, it is better to align EBSD
% with SEM.  This might result in too large matrices due to high resolution
% DIC data, so at this stage, the code needs to be written seperately for
% trace analysis.
%
% The purpose of the code will be mainly, align EBSD to SEM to find grain
% boundaries, therefore to analyze different grains for their traces, or
% neighboring grains traces.  We might need to extend this direction of
% projecting data, but this can be done later.
%
% chenzhe, 2017-06-07.  review code, update some functions.  pay attention
% to FCC.
% note, required functions in chenFunctions
%
% chenzhe, 2017-08-31. change the input format of grain file from .csv to
% .txt, so it is of more general use.
%
% prerequisite: may need EBSD_tool_convert_grain_file_txt_to_CSV()
% ---------- level 1----------
% load_settings()
% grain_file_read()
% find_variable_column_from_grain_file_header()
% % construct_neighbor_structure()
% % construct_misorientation_structure()
% % ----> calculate_misorientation_hcp()
% % --------> euler_to_transformation()
% % --------> derotation()
% find_boundary_from_ID_matrix()
% interp_data()
% generate_grain_avg_data()



clear;
addChenFunction;

[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\Mg4Al_m1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
sampleName = [];    % such as 'Ti7Al_#B6'
cpEBSD = [];    % control points on EBSD image (unit um !!!!)
cpSEM = [];     % control points on SEM image (unit pixel)
sampleMaterial = [];  % such as 'Ti' or 'Mg'
stressTensor = [];
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% data files
[EBSDfileName1, EBSDfilePath1] = uigetfile('E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\EBSD Data\Mg4Al_grain_file_type_1.txt','choose the EBSD file (txt format, from type-1 grain file)');
[EBSDfileName2, EBSDfilePath2] = uigetfile('E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\EBSD Data\Mg4Al_grain_file_type_2.txt','choose the EBSD file (txt format, from type-2 grain file)');
[strainFileName, strainFilePath] = uigetfile('E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\SEM Data\stitched DIC\DIC_merged_3.mat','choose one of the strain file (mat format) for aligning purpose');

% This defines the overlay relationship, ebsdpoint(x,y) * tMatrix = sempoint(x,y)
tform = make_average_transform('projective',cpEBSD,cpSEM);

tMatrix = tform.tdata.T;
tInvMatrix = tform.tdata.Tinv;

saveDataPath = [uigetdir('E:\zhec umich Drive\2021-10-01 MgAl insitu SEM-DIC\Analysis','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
try
    save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],'-append');
catch
    save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
end
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
    
% [temp disable] construct grain neighbor structure S.g1 = [1;2;3], S.g2{ind_in_g1} = [2;3;...]
% neighborStruct = construct_neighbor_structure(EBSDfilePath2,EBSDfileName2);

% [temp disable] Note this is currently only for HCP!
% misorientationStruct.g1 = [1;2;3], misorientationStruct.g2{i} = [2;3;...],
% misorientationStruct.misorientation{i}=[5d;75d;...]
% misorientationStruct = construct_misorientation_structure(neighborStruct, gPhi1, gPhi, gPhi2);


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

%% align euler angle coordinates to SEM, so later on, phiSys can be set to [0 0 0]. Do this before interp so there are few points to rotate
phiSys = [-90, 180, 0];
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); % align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
[q0,q1,q2,q3,phi1,phi,phi2] = regulate_euler_quat(phi1,phi,phi2);   % regulate the angles
[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));
[gQ0,gQ1,gQ2,gQ3,gPhi1,gPhi,gPhi2] = regulate_euler_quat(gPhi1,gPhi,gPhi2);   % regulate the angles
eulerAligned = 1;
save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],'eulerAligned','-append');  % record if eulerAligned.

%% Align EBSD to SEM
% forwarded coordinates 
% [x_EBSD_fwd, y_EBSD_fwd] = tformfwd(tform,x,y);

ID = interp_data(x,y,ID,X,Y,tform,'interp','nearest');
ID(isnan(ID))=0;
phi1 = interp_data(x,y,phi1,X,Y,tform,'interp','nearest');
phi = interp_data(x,y,phi,X,Y,tform,'interp','nearest');
phi2 = interp_data(x,y,phi2,X,Y,tform,'interp','nearest');
edge = interp_data(x,y,edge,X,Y,tform,'interp','nearest');

[gCenterX, gCenterY] = tformfwd(tform,[gCenterX, gCenterY]);

% This part need to find grain average strains.  This should be done after SEM and EBSD is aligned.
% If > 20% data points are valid in this grain, then there will be an avg value.  Otherwise, the value is NaN. ----------------- can modify. 
[gExx,~] = generate_grain_avg_data(ID,gID,exx, 0.2, sigma);

%% special case, if the EBSD data was processed to have gbAdjusted, to align with SEM-DIC, and want to read the adjusted EBSD data
if strcmpi(fileSetting, 'WE43_T6_C1_reducedAOI_setting.mat')
    load('D:\WE43_T6_C1\Analysis_by_Matlab_reducedAOI\EBSD_reduced_AOI.mat');
end

%% Boundary and Unique Boundary Map
[boundaryTF, ~, ~, ~, ~] = find_boundary_from_ID_matrix(ID);
boundaryTFB = boundaryTF;  % make it thicker
for iThick = 1:4
    boundaryTFB = grow_boundary(boundaryTFB);       % grow boundary thicker
end
[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:));
uniqueBoundaryList(uniqueBoundaryList==0) = [];
distMap = bwdist(boundaryTF);

myplot(exx,boundaryTF);
%% neighborStruct and misorinetatinStruct
for ii = 1:length(gID)
    neighborStruct.g1(ii) = gID(ii);
    nbs = gNeighbors(ii,:);
    nbs(nbs==0)=[];
    neighborStruct.g2{ii} = nbs;
end
misorientationStruct = construct_misorientation_structure(neighborStruct, gPhi1, gPhi, gPhi2);

%% Save the data
disp('saving ...');
save([saveDataPath,sampleName,'_traceAnalysis_WS2_rename.mat']);
save([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'], 'ID','X','Y','x','y','boundaryTF','boundaryTFB','eulerAligned',...
    'distMap','uniqueBoundary','uniqueBoundaryList',...
    'phi1','phi','phi2','q0','q1','q2','q3',...
    'gPhi1','gPhi','gPhi2','gQ0','gQ1','gQ2','gQ3',... 
    'neighborStruct','misorientationStruct',...    
    'gArea','gCenterX','gCenterY','gDiameter','gExx','gID','gNNeighbors',...
    'gNeighbors','gPhi1','gPhi','gPhi2','exx','exy','eyy','sigma',...
    'ebsdStepSize','fileSetting','pathSetting',...
    'sampleName','sampleMaterial','stressTensor');
save([saveDataPath,sampleName,'_organized_data'], 'ID','X','Y','x','y','boundaryTF','eulerAligned',...
    'phi1','phi','phi2',...
    'gPhi1','gPhi','gPhi2',... 
    'gArea','gCenterX','gCenterY','gDiameter','gID','gNNeighbors',...
    'gNeighbors');
disp('saving finished');
