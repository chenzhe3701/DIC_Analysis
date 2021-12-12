% chenzhe, 2018-03-30
%
% After converting clustered and manually cleaned twin results, I find that
% some clusters are too small to check.  Maybe it's not good to treat them
% as either twin or nontwin.
% Therefore, in this code, manually provide some stop_grain_clusters to
% check their average pixel size.
% This can be used as thresholds later.
%
% result: shows that maybe threshold = 1000 is a decent choice.

%%
% data required: 'cluster_to_twin_result', load clusterNumMap & _cleaned

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');


% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'ID');

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;



%%
egcs = [201371, 209051, 301372, 303341, 306751, 312551, 400181, 401371, 401373, 405201, 409671, 409672, 503791, 515461];
    
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    cNumMaps{iE} = clusterNumMap;
    cNumMaps_cleaned{iE} = clusterNumMapCleaned;
    struCell{iE} = stru;
end

sz = [];
for ii = 1:length(egcs)

    egc = egcs(ii);
    iC_target = mod(egc,10);
    egc = (egc-iC_target)/10;
    ID_target = mod(egc,10000);
    iE_target = (egc-ID_target)/10000;
    
    ind = (ID==ID_target);
    sz = [sz,sum(cNumMaps{iE_target}(ind) == iC_target)];   


end

figure;
histogram(sz);

% summary: [84, 224, 98, 202, 424, 330, 438, 147, 140, 754, 1262, 1852, 80,
% 148].
%
% So, maybe choose sz = 1000.