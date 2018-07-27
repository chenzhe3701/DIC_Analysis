

% chenzhe, 2018-02-26 based on modified/cleaned cluster result, put cluster
% images into twin and notwin folder. 
% Can use as ground truth, but in fact, it still contain mis-identified
% clusters, because visualization can miss a lot of things ...

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

allImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose/make a parent path for all twin & notwin images.'),'\'];
mkdir(allImgPath,'twin');
mkdir(allImgPath,'notwin');
trainImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose/make a parent path for training twin images.'),'\'];

useCleanedMap = 1;
%% select iE to analyze
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
img_size = 227; % 227 for alexnet, 224 for vgg, googlenet
for iE = iE_start:iE_stop
    % load data for this iE
    warning('off','all');
    %     iE = 5;

    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result]);
    
    if useCleanedMap
        clusterNumMap = clusterNumMapCleaned;
    end
    
    struCell{iE} = stru;
    hWaitbar = waitbar(0,'Converting clustered data into image...');
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
        
        ratio = sz/img_size;
        img = img(ratio:ratio:end,ratio:ratio:end);
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            imgLocal = uint8((img==cNum)*255);
            imgLocal = cat(3,imgLocal,imgLocal,imgLocal);
            
            imgName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
            imwrite(imgLocal, fullfile(allImgPath,imgName));   % to parent folder
            % then to each labeld folder
            if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>=0)
                imwrite(imgLocal, fullfile(allImgPath,'twin',imgName));
            else
                imwrite(imgLocal, fullfile(allImgPath,'notwin',imgName));
            end
        end       
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    try
        close(hWaitbar);
    end
    warning('on','all');
end

%% randomly select nSamples of training image, of twin vs notwin, using current labels  
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    struCell{iE} = stru;
end
mkdir(trainImgPath,'twin');
mkdir(trainImgPath,'notwin');
nSamples = 160;
nTwin = 0;
nNotwin = 0;
while (nTwin<nSamples) || (nNotwin<nSamples)
    iE = randi([iE_start,iE_stop]);
    iS = randi([1,length(stru)]);
    ID_current = gIDwithTrace(iS);
    iCluster = randi([1,length(struCell{iE}(iS).cLabel)]);
    cNum = struCell{iE}(iS).cLabel(iCluster);
    imgName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
    if (struCell{iE}(iS).c2t(iCluster)>nss) && (struCell{iE}(iS).cEnable(iCluster)>=0) && (nTwin<nSamples)
        copyfile(fullfile(allImgPath,imgName), fullfile(trainImgPath,'twin',imgName));
        nTwin = length(dir(fullfile(trainImgPath,'twin')))-2;
    elseif (nNotwin<nSamples)
        copyfile(fullfile(allImgPath,imgName), fullfile(trainImgPath,'notwin',imgName));
        nNotwin = length(dir(fullfile(trainImgPath,'notwin')))-2;
    end
end
