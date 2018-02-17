

% chenzhe, 2018-02-12, make clustered regions into images.
% prepare for cnn work.

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

saveImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] twin image.'),'\'];
%% select iE to analyze

for iE = iE_start:iE_stop
    %% load data for this iE
    warning('off','all');
    %     iE = 5;

    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    load([saveDataPath,name_result_on_the_fly]);
    
    
    %%
    
    hWaitbar = waitbar(0,'Matching cluster with twin system ...');
    for iS =1:length(stru)
        
        ID_current = gIDwithTrace(iS);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        
        [nR,nC] = size(clusterNumMapLocal);
        sz = max(nR,nC);
        img = zeros(sz);
        if nR<nC
            img(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = clusterNumMapLocal;
        else
            img(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = clusterNumMapLocal;
        end
        
        ratio = sz/227;
        img = img(ratio:ratio:end,ratio:ratio:end);
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            imgLocal = uint8((img==cNum)*255);
            imgLocal = cat(3,imgLocal,imgLocal,imgLocal);
            
            imgName = ([num2str(iE*1000000 + ID_current*10 + iCluster),'.tif']);
            
            imwrite(imgLocal, [saveImgPath,imgName]);
        end
        
        
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    try
        close(hWaitbar);
    end
    warning('on','all');
end

%% randomly select nSamples of image, divided them into twin vs notwin by hand later 
trainImgPath = [saveDataPath,'train_img'];
mkdir(trainImgPath);
d=dir(saveImgPath);
d=d(3:end);
nSamples = 4;
ind_pool = randsample(size(d,1),nSamples);
for ii = 1:nSamples
    ind = ind_pool(ii);
    copyfile(fullfile(d(ind).folder,d(ind).name),fullfile(trainImgPath,d(ind).name));
end


