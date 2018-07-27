% for 2018 paper-1
% make barplot for Effective strain vs. distance to grain boundary

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);

dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 6;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);
%% Strain analysis perfroemd to selected strain level.  Can consider save cleaned strain data in the future.  But maybe not needed because it's not difficult to clean and load, although can taken a few seconds.
% Also, the criterion for cleanning strain data might not NEED to be the same? depending on the purpose. 

% Strain data of selected strain level

iE = 5; % target iE

strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
load(strainFile,'exx','exy','eyy','sigma');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------
clear('exy_corrected');
load(strainFile,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.

if exist('exy_corrected','var')&&(1==exy_corrected)
    disp('================= exy already corrected ! ========================');
    exy_corrected = 1;
else
    disp('================= exy being corrected here ! =======================');
    exy = -exy;
    exy_corrected = 1;
end
% remove bad data points
exx(sigma==-1) = nan;
exy(sigma==-1) = nan;
eyy(sigma==-1) = nan;
eEff = calculate_effective_strain(exx,exy,eyy);

% maps for grain boundaries, triple points, and distances to these positions, etc
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'ID','boundaryTF','X','Y','gID','gCenterX','gCenterY','gDiameter','gPhi1','gPhi','gPhi2');  % ID map was pre-interpolated and save.

gEEff = generate_grain_avg_data(ID,gID,eEff,0.01,sigma);
%% reduce resolution for debug ... temp
ratio = 10;
ID_r = ID(1:ratio:end,1:ratio:end);
X_r = X(1:ratio:end,1:ratio:end);
Y_r = Y(1:ratio:end,1:ratio:end);
eEff_r = eEff(1:ratio:end,1:ratio:end);

%%
% which grain it is adjacent to (neighborID)
[boundaryTF, boundaryID, neighborIDonGb, tripleTF, tripleID]...
    = find_boundary_from_ID_matrix(ID_r);

% For each point, calculate its distance to its grain boundary (unit is the same as x,y,gCenterX,gCenterY)
% In addition, find its nearest neighboring grain (i.e., make neighborID complete)
[dToBoundary,dToCenter,dToTriple, rToBoundary,rToCenter,rToTriple, neighborID]...
= find_d_to_boundary_center_triple(X_r,Y_r,ID_r,...
boundaryID,tripleID,neighborIDonGb,gID,...
gCenterX,gCenterY,gDiameter);

phiError = [0 0 0];
if eulerAligned
    phiSys = [0 0 0];
    disp('========== euler already aligned ==============');
else
    phiSys = [-90 180 0];
    disp('======== euler no Aligned, using setting-2 ===========');
end

% [sfBa,sfPr,sfPy,sfPyCA,sfTT,sfBaNum,sfPrNum,sfPyNum,sfPyCANum,sfTTNum]...
%     = calculate_SFs(phi1, phi, phi2, ss, phiSys(1),phiSys(2),phiSys(3), phiError(1),phiError(2),phiError(3), stressTensor, sampleMaterial);

% calculate grain average, then assigne to pixel, [faster]
[gSFBa,gSFPr,gSFPy,gSFPyCA,gSFTT, gSFBaNum,gSFPrNum,gSFPyNum,gSFPyCANum,gSFTTNum, gSF]...
    = calculate_SFs(gPhi1,gPhi,gPhi2, ss, phiSys(1),phiSys(2),phiSys(3), phiError(1),phiError(2),phiError(3), stressTensor, sampleMaterial);

sfBa = assign_field_to_cell(ID_r,gID,gSFBa);
sfPr = assign_field_to_cell(ID_r,gID,gSFPr);
sfPy = assign_field_to_cell(ID_r,gID,gSFPy);
sfPyCA = assign_field_to_cell(ID_r,gID,gSFPyCA);
sfTT = assign_field_to_cell(ID_r,gID,gSFTT);

%% (1) boxplot for eEff vs dToBoundary, pixelwise
data_x = dToBoundary(:)/4096*360; % in micron

bins = [0:10:100,inf];
Labels = zeros(size(data_x,1),1);

legend_text = cell(1,length(bins)-1);
for ii = 1:length(bins)-1
    inds = (data_x(:)>bins(ii));  % they will be corrected again in later loops. 
    Labels(inds) = ii;
    legend_text{ii} =[num2str(bins(ii)),'-',num2str(bins(ii+1))];
end
legend_text{end} = ['>',num2str(bins(end-1))];

% eliminate data that do not belong to any bin
ind = (Labels>0);
Labels = Labels(ind);
data_y = eEff_r(ind);

figure; myboxplot(data_y,Labels,'Whisker',[0.023 0.977]);
set(gca,'ylim',[0,0.15],'fontsize',12,'xticklabel',legend_text,'xticklabelrotation',30);
xlabel('Distance to Grain Boundary, \mum');ylabel('Effective Strain');
print('eEff vs dToGB.tif','-dtiff');

%% (2) boxplot for eEff vs basalSF, pixel wise

data_x = sfBa(:); % in micron

bins = [0:0.05:0.5];
Labels = zeros(size(data_x,1),1);
Labels_c = categorical(size(data_x,1),1);

legend_text = cell(1,length(bins)-1);
for ii = 1:length(bins)-1
    inds = (data_x(:)>bins(ii));  % they will be corrected again in later loops. 
    Labels(inds) = ii;
    legend_text{ii} =[num2str(bins(ii)),'-',num2str(bins(ii+1))];
end

% eliminate data that do not belong to any bin
ind = (Labels>0);
Labels = Labels(ind);
data_y = eEff_r(ind);

figure; myboxplot(data_y,Labels,'Whisker',[0.023 0.977]);
set(gca,'ylim',[0,0.15],'fontsize',12,'xticklabel',legend_text,'xticklabelrotation',30);
xlabel('Basal Schmid Factor');ylabel('Effective Strain');
print('eEff vs sfBa.tif','-dtiff');


%% (3) boxplot for eEff vs twinSF, pixel wise, [ not proportional to sfTT ...]

data_x = sfTT(:); % in micron

bins = [-0.5,0:0.05:0.5];
Labels = zeros(size(data_x,1),1);
Labels_c = categorical(size(data_x,1),1);

legend_text = cell(1,length(bins)-1);
for ii = 1:length(bins)-1
    inds = (data_x(:)>bins(ii));  % they will be corrected again in later loops. 
    Labels(inds) = ii;
    legend_text{ii} =[num2str(bins(ii)),'-',num2str(bins(ii+1))];
end

% eliminate data that do not belong to any bin
ind = (Labels>0);
Labels = Labels(ind);
data_y = eEff_r(ind);

figure; myboxplot(data_y,Labels,'Whisker',[0.023 0.977]);
set(gca,'ylim',[0,0.15],'fontsize',12,'xticklabel',legend_text,'xticklabelrotation',30);
xlabel('Twin Schmid Factor');ylabel('Effective Strain');

%% (4) boxplot for eEff vs sfBa, grain wise

data_x = gSFBa(:); % in micron

bins = [0:0.05:0.5];
Labels = zeros(size(data_x,1),1);
Labels_c = categorical(size(data_x,1),1);

legend_text = cell(1,length(bins)-1);
for ii = 1:length(bins)-1
    inds = (data_x(:)>bins(ii));  % they will be corrected again in later loops. 
    Labels(inds) = ii;
    legend_text{ii} =[num2str(bins(ii)),'-',num2str(bins(ii+1))];
end

% eliminate data that do not belong to any bin
ind = (Labels>0);
Labels = Labels(ind);
data_y = gEEff(ind);

figure; myboxplot(data_y,Labels,'Whisker',[0.023 0.977]);
set(gca,'ylim',[0,0.15],'fontsize',12,'xticklabel',legend_text,'xticklabelrotation',30);
xlabel('Grain Basal Schmid Factor');ylabel('Grain Effective Strain');


%% (5) boxplot for eEff vs sfBa, grain wise

data_x = gSFTT(:); % in micron

bins = [-0.5,0:0.05:0.5];
Labels = zeros(size(data_x,1),1);
Labels_c = categorical(size(data_x,1),1);

legend_text = cell(1,length(bins)-1);
for ii = 1:length(bins)-1
    inds = (data_x(:)>bins(ii));  % they will be corrected again in later loops. 
    Labels(inds) = ii;
    legend_text{ii} =[num2str(bins(ii)),'-',num2str(bins(ii+1))];
end

% eliminate data that do not belong to any bin
ind = (Labels>0);
Labels = Labels(ind);
data_y = gEEff(ind);

figure; myboxplot(data_y,Labels,'Whisker',[0.023 0.977]);
set(gca,'ylim',[0,0.15],'fontsize',12,'xticklabel',legend_text,'xticklabelrotation',30);
xlabel('Grain Twin Schmid Factor');ylabel('Grain Effective Strain');








    
