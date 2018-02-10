

% just look at dissimilarity vs SF.
% chenzhe, 2018-02-09

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);

gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

%% plot dissimilarity vs schmid factor 

for iE = iE_start:iE_stop
    
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    load([saveDataPath,name_result_on_the_fly]);
    
    dissimilarity = [];
    schmidfactor = [];
    for iS=1:length(stru)
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            dissimilarity  = [dissimilarity, pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain)];
            schmidfactor = [schmidfactor, stru(iS).tSF];
            
        end
    end
    figure;plot(dissimilarity,schmidfactor,'.');
    title(['Strain level: ',num2str(STOP{iE+B})]);
    xlabel('Dissimilarity');
    ylabel('Schmid Factor');
    
end
