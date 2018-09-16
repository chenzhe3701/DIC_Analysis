
% relabel twin, separated from previous trace_analysis_2().
% chenzhe, 2018-02-05
% chenzhe, 2018-02-26, based on 2_label_twins_post().

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
end
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

twinImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose the parent path for all twin image.'),'\'];
[fileNet,pathNet] = uigetfile([saveDataPath,'\trainedAlexNet.mat'],'select the trained cnn net');
net = load([pathNet,fileNet],'netTransfer');
net = net.netTransfer;
kk = find(cellfun(@(x) strcmp(x,'twin'),net.Layers(end).ClassNames));
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

useCleanedMap = 0;      % in this code, cleaned or not is only related to the map drawn.  No training, or classification is used in this code.
% But, do pay attention to the image folder you selected.

%% select iE to analyze
for iE = iE_start:iE_stop
    %     iE = 3;
    
    % =========== match cluster with twin system again, if need to change parameter ============================
    % name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];
    % load([saveDataPath,name_result_on_the_fly]);
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    if useCleanedMap
        clusterNumMap = clusterNumMapCleaned;
    end
    % Create a few maps to record the criterion.
    cnnTwinMap = zeros(size(exx));
    
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
        
        cnnTwinMapLocal = zeros(size(ID_local));
        
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);

        stru(iS).tProb = zeros(nCluster,1);
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            indClusterLocal = (clusterNumMapLocal==cNum);
            
            % cnn identify twin. Note that for my processed data, if there is no image, that means the cluster is so small/noisy, so no pixels left for image.     
            imgName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
            try
                I = imread([twinImgPath,imgName]);
            catch
                I=uint8(zeros(227,227,3));
            end
            predictedLabels = classify(net,I);
            predictedProb = predict(net,I);
            
            stru(iS).tProb(iCluster) = predictedProb(kk);
            cnnTwinMapLocal(indClusterLocal) = predictedProb(kk);
        end
        
        % copy identified twin system number to twinMap
        cnnTwinMap(indR_min:indR_max, indC_min:indC_max) = cnnTwinMap(indR_min:indR_max, indC_min:indC_max) + cnnTwinMapLocal;
        
        waitbar(iS/length(stru), hWaitbar);
    end

    
    try
        close(hWaitbar);
    catch
    end
    
    disp('append cnn result to stru');
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    save([saveDataPath,fName_c2t_result],'stru','cnnTwinMap','-append')
end


