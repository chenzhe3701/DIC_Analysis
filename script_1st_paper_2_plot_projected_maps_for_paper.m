
% Seems that this compares the combinations of affine/projective and
% averaged affine/projective to align SEM and EBSD
%
% chenzhe, comment on 2018-01-06.  The code was written on 2017-12-07

clear;
addChenFunction;

[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
sampleName = [];    % such as 'Ti7Al_#B6'
cpEBSD = [];    % control points on EBSD image (unit um !!!!)
cpSEM = [];     % control points on SEM image (unit pixel)
sampleMaterial = [];  % such as 'Ti' or 'Mg'
stressTensor = [];
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% data files
[EBSDfileName1, EBSDfilePath1] = uigetfile('.txt','choose the EBSD file (txt format, from type-1 grain file)');
[EBSDfileName2, EBSDfilePath2] = uigetfile([EBSDfilePath1,'.txt'],'choose the EBSD file (txt format, from type-2 grain file)');
[strainFileName, strainFilePath] = uigetfile([EBSDfilePath1,'.mat'],'choose one of the strain file (mat format) for aligning purpose');

% This defines the overlay relationship, ebsdpoint(x,y) * tMatrix = sempoint(x,y)
% tform = maketform('projective',cpEBSD(1:4,:),cpSEM(1:4,:));
tform = make_average_transform('projective',cpEBSD,cpSEM);

tMatrix = tform.tdata.T;
tInvMatrix = tform.tdata.Tinv;

%%
% Read EBSD grain file. '.txt' format for now.
[EBSDdata1,EBSDheader1] = grain_file_read([EBSDfilePath1, EBSDfileName1]);

columnIndex1 = find_variable_column_from_grain_file_header(EBSDheader1,...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});

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
ID = reshape(EBSDdata1(:,columnIndex1(1)),mResize,nResize)';    ID_0 = ID;
edge = reshape(EBSDdata1(:,columnIndex1(7)),mResize,nResize)';

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
clear EBSDdata1 EBSDdata2 strainData;

%% (1) raw Ebsd-gb plot on original position, i.e., start with (0,0)
cpEBSD;
cpSEM;
fontsize = 24;

tform = make_average_transform('projective',cpEBSD,cpSEM);
boundaryTF_old = find_boundary_from_ID_matrix(ID_0);

myplot(x,y,boundaryTF_old);     % raw ebsd, in um
set(gca,'xlim',[-500,6500],'ylim',[-100,3500])
xlabel('x, micron');ylabel('y, micron');
set(gca,'fontsize',fontsize,'linewidth',2)
%% (2) Ebsd-gb plot on position after aligned with SEM-DIC, i.e., starting from (negative)
[x_EBSD_fwd, y_EBSD_fwd] = tformfwd(tform,x,y);
myplot(x_EBSD_fwd/4096*360, y_EBSD_fwd/4096*360, boundaryTF_old)  % position transformed
set(gca,'xlim',[-500,6500],'ylim',[-100,3500])
xlabel('x, micron');ylabel('y, micron');
set(gca,'fontsize',fontsize,'linewidth',2)

%% (3) Similar to (1), but in black & white
figure;
ind = find(boundaryTF_old(:)==1);
plot3(x(ind), y(ind), boundaryTF_old(ind),'.','markersize',4,'color',[.1 .1 .1])
colorbar;axis equal;
set(gca,'ydir','reverse','xlim',[-500,6500],'ylim',[-100,3500])
xlabel('x, micron');ylabel('y, micron');
set(gca,'fontsize',fontsize,'linewidth',2)
view(0,90)


%% (4) Similar to (2), but Ebsd-gb in black, and overlay with a strain map
skip = 10;
myplot(X(1:skip:end,1:skip:end)/4096*360,Y(1:skip:end,1:skip:end)/4096*360,exx(1:skip:end,1:skip:end));
grid off
hold on;
boundary_surf = ones(size(boundaryTF_old)).*boundaryTF_old*100;
ind = find(boundaryTF_old(:)==1);
plot3(x_EBSD_fwd(ind)/4096*360, y_EBSD_fwd(ind)/4096*360, boundaryTF_old(ind),'.','markersize',4,'color',[.1 .1 .1])
set(gca,'xlim',[-500,6500],'ylim',[-100,3500])
xlabel('x, micron');ylabel('y, micron');
set(gca,'fontsize',fontsize,'linewidth',2)