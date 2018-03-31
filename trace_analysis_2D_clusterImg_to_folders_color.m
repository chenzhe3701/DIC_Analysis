

% chenzhe, 2018-02-26 based on modified/cleaned cluster result, put cluster
% images into twin and notwin folder. 
% Can use as ground truth, but in fact, it still contain mis-identified
% clusters, because visualization can miss a lot of things ...

clear;
addChenFunction;
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

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
    strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile)
    clear('exy_corrected');
    load(strainFile,'exx','exy','eyy','sigma','exy_corrected');     % if 'exy_corrected' does not exist, this does not give error, rather, just warning. % ----------------------------------------------------------------------------------
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
    qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
    qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
    qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
    ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
    exx(ind_outlier) = nan;
    exy(ind_outlier) = nan;
    eyy(ind_outlier) = nan;
    
    % scale the value, prepare for image. !!!
    exx = mat_to_image(exx, [-0.14, 0.07], 'index');
    exy = mat_to_image(exy, [-0.07, 0.07], 'index');
    eyy = mat_to_image(eyy, [-0.07, 0.14], 'index');

    % load data for this iE
    warning('off','all');
    %     iE = 5;

    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','clusterNumMapCleaned');
    
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
        
        exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
        exy_local = exy(indR_min:indR_max, indC_min:indC_max);
        eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        exx_local(ID_local~=ID_current) = 0;
        exy_local(ID_local~=ID_current) = 0;
        eyy_local(ID_local~=ID_current) = 0;
        
        
        [nR,nC] = size(clusterNumMapLocal);
        sz = max(nR,nC);
        img_local_cNum = zeros(sz);
        img_exx_local = zeros(sz);
        img_exy_local = zeros(sz);
        img_eyy_local = zeros(sz);
        if nR<nC
            img_local_cNum(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = clusterNumMapLocal;
            img_exx_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = exx_local;
            img_exy_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = exy_local;
            img_eyy_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = eyy_local;
        else
            img_local_cNum(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = clusterNumMapLocal;
            img_exx_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = exx_local;
            img_exy_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = exy_local;
            img_eyy_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = eyy_local;
        end
        
        ratio = sz/img_size;
        img_local_cNum = img_local_cNum(ratio:ratio:end,ratio:ratio:end);
        img_exx_local = img_exx_local(ratio:ratio:end,ratio:ratio:end);
        img_exy_local = img_exy_local(ratio:ratio:end,ratio:ratio:end);
        img_eyy_local = img_eyy_local(ratio:ratio:end,ratio:ratio:end);
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            % if cVolCleaned==0, move on to the next cluster. The following 3 lines is east to comment out if necessary.
            if (isfield(stru,'cVolCleaned'))&&(stru(iS).cVolCleaned(iCluster)<=0)
                continue;
            end
        
            cNum = stru(iS).cLabel(iCluster);
            % imgLocal = uint8((img_local_cNum==cNum)*255);
            exx_R_Local = uint8(  img_exx_local.*(img_local_cNum==cNum) *255);
            exy_G_Local = uint8(  img_exy_local.*(img_local_cNum==cNum) *255);
            eyy_B_Local = uint8(  img_eyy_local.*(img_local_cNum==cNum) *255);
            
            % imgLocal = cat(3,imgLocal,imgLocal,imgLocal);
            imgLocal = cat(3,exx_R_Local,exy_G_Local,eyy_B_Local);
            
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
nSamples = 120;
nTwin = 0;
nNotwin = 0;
while (nTwin<nSamples) || (nNotwin<nSamples)
    iE = randi([iE_start,iE_stop]);
    iS = randi([1,length(stru)]);
    ID_current = gIDwithTrace(iS);
    for iCluster = (struCell{iE}(iS).cLabel)'
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
end
