
% copy the field stru(ii).tProb, and the map cnnTwinMap,
% from 'cnnTwinResult.mat' to 'cluster_to_twin_result.mat'
%
% chenzhe, 2018-03-01

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;


%%
for iE = iE_start:iE_stop
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cnnTwinResult.mat'];
    load([saveDataPath,fName_c2t_result],'stru','cnnTwinMap');
    stru1 = stru;
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    for ii=1:length(stru)
        stru(ii).tProb = stru1(ii).tProb;
    end
    
    save([saveDataPath,fName_c2t_result],'stru','cnnTwinMap','-append');
    disp(['append saved iE= ',num2str(iE)]);
    
end