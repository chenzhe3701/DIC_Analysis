% chenzhe, 2018-05-14
% In branch 'Ti7Al_B6'
% Organize this to look at grain boundary alignment

clear;
addChenFunction;
dicPath = uigetdir('E:\Ti7Al_B6_insitu_tension\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('E:\Ti7Al_B6_insitu_tension\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'gID','gNeighbors','gNNeighbors','ID','ID_0','iThick','x','y','X','Y','exx','gExx');

gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 14;
% resReduceRatio = 3;         % to save space, reduce map resolution
% grow_boundary_TF = 0;       % whether to grow boundary to make it thicker

% file name prefixes
f1 = 'Ti7Al_B6_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
% twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

%%