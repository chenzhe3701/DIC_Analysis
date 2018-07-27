
% look at how many clusters we have


clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 1;   % elongation levels to analyze. 0-based.
iE_stop = 6;

nC = [];
for iE=iE_start:iE_stop
    name_result_on_the_fly = [sampleName,'_s',num2str(iE+B),'_cluster_result_on_the_fly.mat'];
    load([saveDataPath,name_result_on_the_fly],'stru');
    
    for iS=1:length(stru)
       nC{iE}(iS) = length(stru(iS).cLabel);
    end
    figure; histogram(nC{iE});

end

