% temporary code with 2 functions:
% 
% (1) original purpose, associate twin variant with grain boundary.
%
% (2) also, maybe can use it to divide clusterMap into trueTwinMap.


clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_2021_09\','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\WE43_T6_C1\Analysis_2021_09\20211001_0212_relabeled_result.mat','select the confirmed labeled results, for struCell');

% This provides a valid struCell
% [twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the results for twin-grain boundary intersection');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
end
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

debugTF = 0;

% [data] strain data. Convert into v7.3 for partial loading
clear strainFile;
for iE = iE_start-1:iE_stop
    strainFileName = [dicPath,'\',f2,STOP{iE+B}];
    disp(strainFileName);
    if ~exist([strainFileName,'_v73.mat'],'file')
        load(strainFileName);
        clear('exy_corrected');
        load(strainFileName,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
        
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
        
        outlier_removed = 1;
        save([strainFileName,'_v73.mat'],'outlier_removed','exy_corrected','-v7.3');
        
        myFile = matfile(strainFileName);
        myFields = who(myFile);
        for ii=1:length(myFields)
            save([strainFileName,'_v73.mat'],myFields{ii},'-append','-v7.3');
        end
    else
        disp('v7.3 file already exist');
    end
    strainFile{iE} = matfile([strainFileName,'_v73.mat']);
end

% (0) load data, using SF threshold values to assign active twin system, and make maps
% matfile that contains (cluster number maps)  
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    c2tFile{iE} = matfile([saveDataPath,fName_c2t_result]);
    clusterNumMapCell{iE} = c2tFile{iE}.clusterNumMapCleaned;
end

% Load 'struCell'.  This script will assign
% 'trueTwinMapCell'='variantMapCleanedCell'
% 'cToGbDistMapCell' is omitted, as will no longer be used in this code)
load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell');   

% load previous twin_gb interaction result, for reference.
% load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile),'struCell');

%%
% Assign each variant to each boundary.
%
% Algorithm:
% (1) crop [X,Y,uniqueGB,trueTwinVariant]_local map.
%
% (2) For each variant, rotate the maps so that the theoretical trace
% direction is horizontal.
%
% (3) For each row, find two intersecting gbs [gbL, gbR].  Starting from
% both ends, iteratively assign each connected segment to gbL or gbR,
% depending on the gbLabel of the closest neighbor pt on that line within
% that variant. (Record it to matrix gb_belonged_to).
% 
% (4) Create a distance ruler to help. If a boundary has >XXX(e.g.,8)
% points within distance XXX(e.g,10), or within 50% of the segment length
% if the horizontal line is < 20 pts long, it has an intersection point to
% the variant.  Record/label the gb as having_intersection.
% 
% (5) Go back, and check each row's data points' two associated gbLabels.

for iE = 2:5
    variantMap = zeros(size(ID));
    for iS = 1:length(struCell{iE}) % 23
        ID_current = struCell{iE}(iS).gID;
        disp(['check iE=',num2str(iE),', ID=',num2str(ID_current)]);
        ind = find(gID==ID_current);
        
        euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        if (1==eulerAligned)
            g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
        else
            g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
        end
        [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
        traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);    % angle x-to-y  
        
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        % (Step-1) Crop local maps
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        X_local = X(indR_min:indR_max, indC_min:indC_max);
        Y_local = Y(indR_min:indR_max, indC_min:indC_max);
        uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max); 
        uniqueBoundary_local((floor(uniqueBoundary_local/10000)~=ID_current)&(mod(uniqueBoundary_local,10000)~=ID_current)) = 0;    % leave only associated with this grain.
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
%         trueTwinMap_local = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);  
%         trueTwinMap_local(ID_local~=ID_current) = 0;
        clusterNumMap_local = clusterNumMapCell{iE}(indR_min:indR_max, indC_min:indC_max);  
        clusterNumMap_local(ID_local~=ID_current) = 0;
        % Find active system, if any, using cTrueTwin/tGb field
        activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        
        grains_variant_map = zeros(size(ID_local));
        % for each cluster, depending on how many variants it has
        for iC = 1:length(struCell{iE}(iS).cLabel)
            if sum(struCell{iE}(iS).cTrueTwin(iC,:))==0
                % not a twinned cluster
                tnMap = zeros(size(clusterNumMap_local));
            elseif sum(struCell{iE}(iS).cTrueTwin(iC,:))==1
                % assign this cluster to a single variant 
                variant_num = find(struCell{iE}(iS).cTrueTwin(iC,:));
                tnMap = zeros(size(clusterNumMap_local));
                tnMap(clusterNumMap_local==iC) = variant_num;
            elseif sum(struCell{iE}(iS).cTrueTwin(iC,:))>1
                clear csl;
                % need to divide
                for iTwin = 1:6
                    if struCell{iE}(iS).cTrueTwin(iC,iTwin)==1
                        % (Step-2) Rotate maps
                        ID_r = imrotate(ID_local, traceDir(iTwin), 'nearest', 'loose');
                        X_r = imrotate(X_local, traceDir(iTwin), 'nearest', 'loose');
                        Y_r = imrotate(Y_local, traceDir(iTwin), 'nearest', 'loose');
%                         uniqueBoundary_r = imrotate(uniqueBoundary_local, traceDir(iTwin), 'bilinear', 'loose');
                        uniqueBoundary_r = imrotate(uniqueBoundary_local, traceDir(iTwin), 'nearest', 'loose');
                        uniqueBoundary_r = imdilate(uniqueBoundary_r,ones(3));
                        clusterNumMap_r = imrotate(clusterNumMap_local, traceDir(iTwin), 'nearest', 'loose');
                        vMap_r = (clusterNumMap_r == iC); % try this to segment clusterNumMap into variantMap(or trueTwinMap)
                        
                        [nr,nc] = size(vMap_r);
                        gbLabelMap = zeros(nr,nc);    % to store assigned gb_label
                        gbNumXY_intersect = [];     % [gbNum, Xpos, Ypos]  ----------------------------------> this could be recorded in struCell.
                        cslMap = zeros(nr,nc);     % a map recording Connected Segment Length (CSL)   ----------> this might be helpful do determine variant number, if starting point is clusterNumMap.
                        gbLR = zeros(nr,2);     % store each rows two possible gbs
                        
                        % (Step-3)
                        for ir = 1:nr
                            if any(vMap_r(ir,:))
                                icL_back = find(uniqueBoundary_r(ir,:),1,'first');
                                gbL = uniqueBoundary_r(ir,icL_back);
                                icR_back = find(uniqueBoundary_r(ir,:),1,'last');
                                gbR = uniqueBoundary_r(ir,icR_back);
                                gbLR(ir,:) = [gbL, gbR];
                                
                                % (instert Step-4) determine if gbL/R can be considered as an intersecting gb. -----------------------------------
                                length_cr = round(min(30, (icR_back - icL_back)/2));
                                num_cr = round(length_cr * 0.7);
                                if sum(vMap_r(ir,icL_back:icL_back+length_cr))>num_cr
                                    gbNumXY_intersect = [gbNumXY_intersect; gbL, X_r(ir,icL_back), Y_r(ir,icL_back)];
                                end
                                if sum(vMap_r(ir,icR_back-length_cr:icR_back))>num_cr
                                    gbNumXY_intersect = [gbNumXY_intersect; gbR, X_r(ir,icR_back), Y_r(ir,icR_back)];
                                end
                                % end of (Step-4). Alternatively, read from manual label for data analysis. ---------------------------------------
                                
                                icL_front = find(vMap_r(ir,:),1,'first');
                                icR_front = find(vMap_r(ir,:),1,'last');
                                % will be false, if either is empty
                                while (icL_front<=icR_front)
                                    csl_length = 0;
                                    if (icL_front-icL_back)<=(icR_back-icR_front)
                                        % search for connected segments from left to right
                                        while(vMap_r(ir,icL_front))
                                            gbLabelMap(ir,icL_front) = 1;  % prepare to assign label
                                            vMap_r(ir,icL_front) = 0; % make element on variant map 0
                                            icL_front = icL_front + 1;   % move pointer forward to the right
                                            csl_length = csl_length + 1;
                                        end
                                        cslMap(gbLabelMap==1) = csl_length;
                                        gbLabelMap(gbLabelMap==1) = gbL;
                                        icL_back = icL_front - 1;    % assign left side back
                                        icL_front = find(vMap_r(ir,:),1,'first'); % search left side front again
                                    else
                                        while (vMap_r(ir,icR_front))
                                            gbLabelMap(ir,icR_front) = 1;
                                            vMap_r(ir,icR_front) = 0;
                                            icR_front = icR_front - 1;
                                            csl_length = csl_length + 1;
                                        end
                                        cslMap(gbLabelMap==1) = csl_length;
                                        gbLabelMap(gbLabelMap==1) = gbR;
                                        icR_back = icR_front + 1;
                                        icR_front = find(vMap_r(ir,:),1,'last');
                                    end
                                end
                            end % end of if any(variant_r(ir,:))
                        end % end of for ir=1:nr
                        
                        % (Step-5) go back to clean.  Sometimes, no intersection is determined ...    
                        cleanTF = 1;
                        if (cleanTF)&&(~isempty(gbNumXY_intersect))
                            gbList = unique(gbNumXY_intersect(:,1));
                            for ir=1:nr
                                tf = ismember(gbLR(ir,:),gbList);
                                if sum(tf)==1
                                    gbOK = gbLR(ir, tf);
                                    gbKO = gbLR(ir, ~tf);
                                    ind = gbLabelMap(ir,:) == gbKO;
                                    gbLabelMap(ir,ind) = gbOK;
                                elseif sum(ismember(gbLR(ir,:),gbList))==0
                                    gbLabelMap(ir,:) = 0;
                                end
                            end
                        end
                        
                        temp = cslMap;
                        % rotate back, need to crop again.
                        temp = imrotate(cslMap,-traceDir(iTwin), 'nearest', 'loose');
                        ID_back = imrotate(ID_r,-traceDir(iTwin), 'nearest', 'loose');
                        ind_back = ismember(ID_back, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
                        
                        % Need to crop a region of the original size.  
                        img1_template = (ID_local==ID_current);
                        img2_signal = (ID_back==ID_current);
                        [yOffSet, xOffSet] = normxcorr2A_register(img1_template, img2_signal, [0 0 0 0], [0 0 0 0], 0);
                        indC_back_min = 1 + xOffSet;
                        indC_back_max = indC_back_min + indC_max - indC_min;
                        indR_back_min = 1 + yOffSet;
                        indR_back_max = indR_back_min + indR_max - indR_min;
                        
                        % Make it one data point wider on each side
%                         indC_back_min = max(1, find(sum(ind_back, 1), 1, 'first')-1);
%                         indC_back_max = min(size(ID_back,2), find(sum(ind_back, 1), 1, 'last')+1);
%                         indR_back_min = max(1, find(sum(ind_back, 2), 1, 'first')-1);
%                         indR_back_max = min(size(ID_back,1), find(sum(ind_back, 2), 1, 'last')+1);
                        
                        csl(:,:,iTwin) = temp(indR_back_min:indR_back_max, indC_back_min:indC_back_max);
                        
                    end % end of if(activeTS(iTwin)==1)
                    
                end % end of for iTwin=1:6
            
                [nr,nc,np] = size(csl);
                tnMap = zeros(nr,nc);   % variant_num_map of this cluster
                for ir=1:nr
                    for ic=1:nc
                        [maxV,temp] = max(csl(ir,ic,:));
                        if maxV>0
                            tnMap(ir,ic) = temp;
                        end
                    end
                end
            
            end % end of elseif sum(struCell{iE}(iS).cTrueTwin(iC,:))>1
            
            % Here need to prevent elements belonging to other clusters having a tnMap value (caused by rotation).   
            tnMap(clusterNumMap_local~=iC) = 0;            
            
            clusters_variant_map{iC} = tnMap;
            grains_variant_map = grains_variant_map + clusters_variant_map{iC};
        end % end of for iC = 1:num_of_clusters
        variantMap(indR_min:indR_max, indC_min:indC_max) = variantMap(indR_min:indR_max, indC_min:indC_max) + grains_variant_map;

    end % end of (for iS=1:end)
    variantMapCell{iE} = variantMap;
end % end of (for iE=2:5) 
% save(fullfile(saveDataPath, [sampleName,'_new_variant_map.mat']),'variantMapCell','-v7.3')  
 
%% try to clean variantMap again 
disp('clean variant map ...');
useParallel = 1;
if useParallel    
    for iE = iE_start:iE_stop
        variantMap = variantMapCell{iE};
        
        npool = 3;
        tN = 1:length(struCell{iE});
        a = npool - mod(length(tN),npool);
        tN = [zeros(1,a),tN];
        tN = reshape(tN,npool,[]);
        ipool = [];
        partMap = [];
        for ii = 1:npool
            ipool{ii} =  tN(ii,:);
            partMap{ii} = zeros(size(ID));
        end
        parfor ii = 1:npool
            for iS =ipool{ii}
                % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
                % iS = find(gIDwithTrace == 296); % for debugging.
                if iS > 0
                    ID_current = struCell{iE}(iS).gID;
                    
                    
                    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
                    indC_min = find(sum(ind_local, 1), 1, 'first');
                    indC_max = find(sum(ind_local, 1), 1, 'last');
                    indR_min = find(sum(ind_local, 2), 1, 'first');
                    indR_max = find(sum(ind_local, 2), 1, 'last');
                    
                    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
                    
                    variantMapLocal = variantMap(indR_min:indR_max, indC_min:indC_max);
                    variantMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
                    
                    variantMapLocal = one_pass_fill(variantMapLocal);
                    
                    partMap{ii}(indR_min:indR_max, indC_min:indC_max) = partMap{ii}(indR_min:indR_max, indC_min:indC_max) + variantMapLocal;
                    disp(['cleaning, ID = ',num2str(ID_current)]);
                end
            end
        end
        
        variantMapCleaned = zeros(size(ID)); % new map of interest
        for ii=1:npool
            variantMapCleaned = variantMapCleaned + partMap{ii};
        end
        
        variantMapCleanedCell{iE} = variantMapCleaned;
    end    
end
save(fullfile(saveDataPath, [sampleName,'_new_variant_map.mat']),'variantMapCleanedCell','-append');


%% Count tVol again, need to load a valid struCell here. ================================================  
SF_th = 0.15;  
for iE = iE_start:iE_stop
    % use the updated trueTwinMapCell.  Use variantMapCleanedCell--------------------   
    trueTwinMapCell{iE} = variantMapCleanedCell{iE};
    trueTwinMapCell{iE}(trueTwinMapCell{iE}>0) = trueTwinMapCell{iE}(trueTwinMapCell{iE}>0) + 18;
    
    ind = trueTwinMapCell{iE}(:)~=0;
    [idic,IA,IC] = unique(ID(ind)*1000+trueTwinMapCell{iE}(ind));   % without using 'ind', it can include point that does not belong to any twin !!! 2018-11-19
    idicCount = histcounts(ID(ind)*1000+trueTwinMapCell{iE}(ind),[idic(:);idic(end)+1]);
    
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        disp(['count tVol, iE=',num2str(iE),', ID=',num2str(ID_current)]);
        struCell{iE}(iS).tVol = zeros(1, length(struCell{iE}(iS).tLabel));  % reset and count again.
        
        for iTwin = 1:length(struCell{iE}(iS).tSF)
            if struCell{iE}(iS).tSF(iTwin) < SF_th 
                struCell{iE}(iS).cTrueTwin(:,iTwin) = 0;
            end
            ind = find(idic==ID_current*1000+struCell{iE}(iS).tLabel(iTwin));
            if ~isempty(ind) % && sum(struCell{iE}(iS).cTrueTwin(:,iTwin))>0
                struCell{iE}(iS).tVol(iTwin) = idicCount(ind);
            end
        end
    end
    
end


% check
struCell{5}(105)
struCell{5}(106)
struCell{5}(107)

save(fullfile(saveDataPath, [sampleName,'_new_variant_map.mat']),'trueTwinMapCell','struCell','-append');
copyfile(fullfile(saveDataPath, [sampleName,'_new_variant_map.mat']),fullfile(saveDataPath, [sampleName,'_final_variant_map.mat']));






%% previous thought, but may not be good
% E.g., v_1, accounts for 40% of grain size.
% 20% associated with gb A, 80% with gb B, 0% with gb C.
% So, 'membership' = 20%, 80%.
% Area fraction distributed to each gb = 8%, 32% --> this can be used in
% analysis of twin-grain boundary interaction, where 8% and 32% can be
% normalized by gb_length.
% This may also be used as an alternative method to determine twin-gb
% intersection? As this can be done after twin variant identification.
