% chenzhe, 2018-03-03
%
% threshold clusterSize
% for cluster/clusterCleaned with size < threshold, 
% (1) change cluster label/clusterCleaned label to  negative.
% (2) if clusterCleaned < threshold, make 'cVol' field as 0.

%%
% data required: 'cluster_to_twin_result', load clusterNumMap (& _cleaned),
% 'stru' to edit the 'cVol' field.

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

threshold = 1000;
%% data
for iE = iE_start:iE_stop
    fName_source = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    copyfile([saveDataPath,fName_source],[saveDataPath,fName_c2t_result],'f');
    
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    
    clusterNumMap = abs(clusterNumMap);
    clusterNumMapCleaned = abs(clusterNumMapCleaned);
    for iS = 1:length(stru)
        
        % initialize
        stru(iS).cVol = zeros(size(stru(iS).cLabel));   
        stru(iS).cVolCleaned = zeros(size(stru(iS).cLabel)); 
        
        for iCluster = 1:size(stru(iS).cLabel)
            cNum = stru(iS).cLabel(iCluster);
            
            ind = (ID==stru(iS).gID)&(abs(clusterNumMap)==cNum);
            sz = sum(ind(:));
            if sz > threshold
                stru(iS).cVol(iCluster) = sz;
            else
                stru(iS).cVol(iCluster) = -sz;
                clusterNumMap(ind) = -cNum;
                disp(['small cluster detected: ID=',num2str(stru(iS).gID),', cluster=',num2str(cNum),', size=',num2str(sz)]);
            end
            
            ind = (ID==stru(iS).gID)&(abs(clusterNumMapCleaned)==cNum);
            sz = sum(ind(:));
            if sz > threshold
                stru(iS).cVolCleaned(iCluster) = sz;
            else
                stru(iS).cVolCleaned(iCluster) = -sz;
                clusterNumMapCleaned(ind) = -cNum;
                disp(['small cluster_cleaned detected: ID=',num2str(stru(iS).gID),', cluster=',num2str(cNum),', size=',num2str(sz)]);
            end
        end
        if (rem(iS,100)==0)
            disp(['finished ',num2str(iS),' grains']);
        end
    end
    save([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned','-append');
end

