% script to get and save local strain data for a specific grain
% chenzhe, 20180308

addChenFunction;

grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
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
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx');



% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';


% The strain data part
clear data;
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    
    
    strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
    load(strainFile,'exx','exy','eyy','sigma');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------
    clear('exy_corrected');
    load(strainFile,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
    
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
    
    %%%% grain
    
    ID_current=ID_target;
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
    exy_local = exy(indR_min:indR_max, indC_min:indC_max);
    eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    sigma_local = sigma(indR_min:indR_max, indC_min:indC_max);
    
    clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
    exx_local(ID_local~=ID_current) = nan;  % The nans are used for making images for the presentation. 0 is OK in processing, because
    exy_local(ID_local~=ID_current) = nan;
    eyy_local(ID_local~=ID_current) = nan;
    sigma_local(ID_local~=ID_current) = nan;
    
    s.iE = iE;
    s.indR_min = indR_min;
    s.indR_max = indR_max;
    s.indC_min = indC_min;
    s.indC_max = indC_max;
    s.boundaryTF_local = boundaryTF_local;
    s.x_local = x_local;
    s.y_local = y_local;
    s.ID_local = ID_local;
    s.ID_current = ID_current;
    s.exx_local = exx_local;
    s.exy_local = exy_local;
    s.eyy_local = eyy_local;
    s.sigma_local = sigma_local;
    s.clusterNumMapLocal = clusterNumMapLocal;
    
    % cleaned
    clusterNumMapCleanedLocal = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapCleanedLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
    s.clusterNumMapCleanedLocal = clusterNumMapCleanedLocal;

    data(iE) = s;

end

fname = [f1,'_all_grain_',num2str(ID_current),'_local_map.mat'];
save(fullfile(grainDataPath,fname),'data');