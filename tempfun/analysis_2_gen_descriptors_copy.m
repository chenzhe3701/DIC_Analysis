% chenzhe 2016-5-23
% analysis_2, Calculate lots of descriptors:
% (1) generate dToGB, boundaryTF, ...
% (2) check: plot strain distribution vs. dToGB

%% load data
addChenFunction;
if ~exist('mainFolder','var')
    mainFolder = uigetdir('','choose the folder where you saved the aligned data');
end
% mainFolder ='H:\Ti7Al #B6 data aligned 2k res';
EBSD_data = [mainFolder,'\EBSD_Data_Aligned_Layer_1'];
load(EBSD_data);

Strain_data = [mainFolder, '\',uigetfile(mainFolder,'choose the strain data')];
% Strain_data = [mainFolder,'\Ti7Al_#B6_new_stop_015'];
load(Strain_data);
e0 =  sqrt(2/3*(exx.^2+eyy.^2+2*exy.^2));
[gE0,~] = generate_grain_avg_data(ID,gID,e0, 0.2, sigma);
save(Strain_data, 'e0', 'gE0', '-append');

%% (1.1) generate uniqueGB and dToGB(manhattan)
uniqueGB = construct_unique_boundary_ID(ID);
boundaryTF = double(logical(uniqueGB));
dToGB_Man = city_block(uniqueGB,-1);   % Manhattan distance
dToBoundary_Man = dToGB_Man*(y(2)-y(1));        % in um

% (1.2) regulate euler angle, and make quarternion description
[q0,q1,q2,q3,phi1New,phiNew,phi2New] = regulate_euler_quat(phi1,phi,phi2);
[gQ0,gQ1,gQ2,gQ3,gPhi1New,gPhiNew,gPhi2New] = regulate_euler_quat(gPhi1,gPhi,gPhi2);
for ii=1:length(gPhi1)
    g = euler_to_transformation([gPhi1(ii),gPhi(ii),gPhi2(ii)], [-90, 180, 0], [0 0 0]);
    [eulerd,dcm,quat,rod,V,thetad] = consider_hcp_symmetry(g);
    rodV(ii,:) = rod;
    Quaternion(ii,:) = quat;
end

% (1.3) calculate Schmid Factors
[ss, c_a] = define_SS('Ti','twin');       % Ti or Mg --------------------------------------------------
ss = hex_to_cart_ss(ss, c_a);
stressTensor = [-1 0 0; 0 0 0; 0 0 0];       % stress, tensile or compressive ----------------------------
[gSFBa,gSFPr,gSFPy,gSFPyCA,gSFTT, gSFBaNum,gSFPrNum,gSFPyNum,gSFPyCANum,gSFTTNum, gSF]...
    = calculate_SFs(gPhi1,gPhi,gPhi2, ss, -90,180,0, 0,0,0, stressTensor, 'HCP-Mg');

sfBa = assign_field_to_cell(ID,gID,gSFBa);
sfPr = assign_field_to_cell(ID,gID,gSFPr);
sfPy = assign_field_to_cell(ID,gID,gSFPy);
sfPyCA = assign_field_to_cell(ID,gID,gSFPyCA);
sfTT = assign_field_to_cell(ID,gID,gSFTT);

% (1.4) calculate distance to grainboundaries and triplepoints
[boundaryTF, boundaryID, neighborIDonGb, tripleTF, tripleID] = find_boundary_from_ID_matrix(ID);

% For each point, calculate its distance to its grain boundary (unit is the same as x,y,gCenterX,gCenterY)
% In addition, find its nearest neighboring grain (i.e., make neighborID complete)
[dToBoundary,dToCenter,dToTriple,rToBoundary,rToCenter,rToTriple, neighborID]...
    = find_d_to_boundary_center_triple(x,y,ID,boundaryID,tripleID,neighborIDonGb,gID,gCenterX,gCenterY,gDiameter);

save(EBSD_data,'uniqueGB','boundaryTF','dToGB_Man','dToBoundary_Man',...
    'q0','q1','q2','q3','phi1New','phiNew','phi2New',...
    'gQ0','gQ1','gQ2','gQ3','gPhi1New','gPhiNew','gPhi2New','rodV','Quaternion',...
    'gSFBa','gSFPr','gSFPy','gSFPyCA','gSFTT','gSFBaNum','gSFPrNum','gSFPyNum','gSFPyCANum','gSFTTNum','gSF',...
    'sfBa','sfPr','sfPy','sfPyCA','sfTT',...
    'tripleTF','dToBoundary','dToCenter','dToTriple','rToBoundary','rToCenter','rToTriple','neighborID','-append');

%% (2) check and plot relationship between strain and dToGB
qts = quantile(e0(:), 0:0.01:1);
ind = (e0>qts(1))&(e0<qts(101))&(~isnan(e0))&(dToBoundary>0.1)&(dToBoundary<30);
eData = e0(ind);
dData = dToBoundary(ind);
d2Data = cell(size(dData));

xPlot = [];yPlot=[];
step = 2;
xVal = [0:step:10,inf];   % last value should be larger than extreme
for ii = 1:length(xVal)-1
    x1 = xVal(ii);
    x2 = xVal(ii+1);
    inds = (dData>=x1)&(dData<x2);
    dData(inds) = x1;
    numDataPoints(ii) = sum(inds);
    for jj=1:length(inds)
        if inds(jj)>0
            d2Data{jj}=[num2str(x1),' - ',num2str(x2)];
        end
    end
    nameOrder{ii} = [num2str(x1),' - ',num2str(x2)];
end
figure; myboxplot(eData,d2Data,'whisker',[0.03,0.97],'grouporder',nameOrder);
set(gca,'ylim',[0,0.15]);
xlabel('Distance to Grain Boundary, micron');ylabel('Effective Strain');
set(gca,'fontsize',16)
%% (3) check and plot relationship between strain and Schmid factor, such as basal, point wise
qts = quantile(e0(:), 0:0.01:1);
ind = (e0>qts(1))&(e0<qts(101))&(~isnan(e0))&(dToBoundary>0.1)&(dToBoundary<30);
eData = e0(ind);
dData = sfBa(ind);      % ----------------- this is the data vs strain. ------------
d2Data = cell(size(dData));

xPlot = [];yPlot=[];
step = 0.1;
xVal = [0:step:0.5];   % last value should be larger than extreme
clear nameOrder;
for ii = 1:length(xVal)-1
    x1 = xVal(ii);
    x2 = xVal(ii+1);
    inds = (dData>=x1)&(dData<x2);
    dData(inds) = x1;
    numDataPoints(ii) = sum(inds);
    for jj=1:length(inds)
        if inds(jj)>0
            d2Data{jj}=[num2str(x1),' - ',num2str(x2)];
        end
    end
    nameOrder{ii} = [num2str(x1),' - ',num2str(x2)];
end
figure; myboxplot(eData,d2Data,'whisker',[0.03,0.97],'grouporder',nameOrder);
set(gca,'ylim',[0,0.15]);
xlabel('Basal Schmid Factor');ylabel('Effective Strain');
set(gca,'fontsize',16)
%% (4) check and plot relationship between gE0 and gDiameter, grain wise

ind = ~isnan(gE0);
eData = gE0(ind);
dData = gDiameter(ind);      % ----------------- this is the data vs strain. ------------
d2Data = cell(size(dData));

xPlot = [];yPlot=[];
step = 5;
xVal = [0:step:25,100];   % last value should be larger than extreme
clear nameOrder;
for ii = 1:length(xVal)-1
    x1 = xVal(ii);
    x2 = xVal(ii+1);
    inds = (dData>=x1)&(dData<x2);
    dData(inds) = x1;
    numDataPoints(ii) = sum(inds);
    for jj=1:length(inds)
        if inds(jj)>0
            d2Data{jj}=[num2str(x1),'-',num2str(x2)];
        end
    end
    nameOrder{ii} = [num2str(x1),'-',num2str(x2)];
end
figure; myboxplot(eData,d2Data,'grouporder',nameOrder);
% set(gca,'ylim',[0,0.12]);
xlabel('Grain Diameter, micron');ylabel('Grain Effective Strain');
set(gca,'fontsize',16);

%% (5) check and plot relationship between gE0 and gSFBa, grain wise

ind = ~isnan(gE0);
eData = gE0(ind);
dData = gSFBa(ind);      % ----------------- this is the data vs strain. ------------
d2Data = cell(size(dData));

xPlot = [];yPlot=[];
step = 0.1;
xVal = [0:step:0.5];   % last value should be larger than extreme
clear nameOrder;
for ii = 1:length(xVal)-1
    x1 = xVal(ii);
    x2 = xVal(ii+1);
    inds = (dData>=x1)&(dData<x2);
    dData(inds) = x1;
    numDataPoints(ii) = sum(inds);
    for jj=1:length(inds)
        if inds(jj)>0
            d2Data{jj}=[num2str(x1),' - ',num2str(x2)];
        end
    end
    nameOrder{ii} = [num2str(x1),' - ',num2str(x2)];
end
figure; myboxplot(eData,d2Data,'grouporder',nameOrder);
% set(gca,'ylim',[0,0.12]);
xlabel('Grain Basal Schmid Factor');ylabel('Grain Effective Strain');
set(gca,'fontsize',12);

%% (6) check and plot relationship between gE0 and gSFTT, for gSFBa<0.2 grain wise

ind = ~isnan(gE0)&(gSFBa<0.2);
eData = gE0(ind);
dData = gSFTT(ind);      % ----------------- this is the data vs strain. ------------
d2Data = cell(size(dData));

xPlot = [];yPlot=[];
step = 0.1;
xVal = [-0.5:step:0.5];   % last value should be larger than extreme
clear nameOrder;
for ii = 1:length(xVal)-1
    x1 = xVal(ii);
    x2 = xVal(ii+1);
    inds = (dData>=x1)&(dData<x2);
    dData(inds) = x1;
    numDataPoints(ii) = sum(inds);
    for jj=1:length(inds)
        if inds(jj)>0
            d2Data{jj}=[num2str(x1),' - ',num2str(x2)];
        end
    end
    nameOrder{ii} = [num2str(x1),' - ',num2str(x2)];
end
figure; myboxplot(eData,d2Data,'grouporder',nameOrder);
% set(gca,'ylim',[0,0.12]);
xlabel('Grain Basal Schmid Factor');ylabel('Grain Effective Strain');
set(gca,'fontsize',12);


%% () continue, if we need traceStructure file
% [traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
% traceStruct = load([traceStructPath,'\',traceStructFile]);
% traceStruct = traceStruct.traceStruct;
% 
% 
% 


