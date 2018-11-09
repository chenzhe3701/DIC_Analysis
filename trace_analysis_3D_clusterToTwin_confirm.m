
% chenzhe, 2018-10-20
%
% After we got an idea of suitable SF

clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[preLabelFile, preLabelPath] = uigetfile('D:\p\m\DIC_Analysis\','select the results where twin identification was based on trace dir and strain');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
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

%% [data] strain data. Convert into v7.3 for partial loading
clear strainFile;
for iE = iE_start:iE_stop
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

%% (0) load data, using SF threshold values to assign active twin system, and make maps
% Load cluster number maps (cleaned).
clusterNumberMapCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
load(fullfile(preLabelPath,preLabelFile),'struCell','twinMapCell','sfMapCell');
% And initialize fields in the 'struCell'
for iE = iE_start:iE_stop
    % initialize/zero related fields
    for iS =1:length(struCell{iE})
        struCell{iE}(iS).cTrueTwin = struCell{iE}(iS).cActiveSS;
        struCell{iE}(iS).tVol = zeros(1, length(struCell{iE}(iS).tLabel));
    end
end

%% First step is to use SF as a threshold, and create/update trueTwinMap, and struCell{iE}.cTrueTwin
trueTwinMapCell = cell(1,length(STOP)-1);
SF_th = 0.15;    % at critical value of 0, already seem to be very helpful.

for iE = iE_start:iE_stop
    % Making a trueTwinMapCell by the SF_th is easy.
    trueTwinMapCell{iE} = twinMapCell{iE};
    trueTwinMapCell{iE}(sfMapCell{iE}<SF_th)=0;
    % But have to modify each field in the struCell
    
    [idic,IA,IC] = unique(ID(:)*1000+trueTwinMapCell{iE}(:));
    idicCount = histcounts(ID(:)*1000+trueTwinMapCell{iE}(:),[idic(:);idic(end)+1]);
    
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        
        for iTwin = 1:length(struCell{iE}(iS).tSF)
            if struCell{iE}(iS).tSF(iTwin) < SF_th
                struCell{iE}(iS).cTrueTwin(:,iTwin) = 0;
            end
            ind = find(idic==ID_current*1000+struCell{iE}(iS).tLabel(iTwin));
            if ~isempty(ind)
                struCell{iE}(iS).tVol(iTwin) = idicCount(ind);
            end
        end
    end
    
end

%% Then we need to plot strain maps, cluster number maps, temporary trueTwinMaps to modify
iE_select = 5;  % select an iE

% plot maps to check: strain map, twin map, cluster number map
% close all;
myplot_twin_id(strainFile{iE_select}.exx,trueTwinMapCell{iE_select},clusterNumberMapCell{iE_select},'tf',boundaryTFB,'x',X,'y',Y);

%% Then [relabel] and [recalculated tVol].
warning('off','all');
for ii = 2:size(tNote,1)
    ID_current = tNote(ii,1);
    iE = tNote(ii,2);
    iC = tNote(ii,3);
    activeSS = tNote(ii,4:9);
    
    iS = find(arrayfun(@(x) x.gID == ID_current,struCell{iE}));
    if ID_current > 0
        
        % (1) Calculate theoretical trace direction.
        ind_euler = find(gID==ID_current);
        euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
        if (1==eulerAligned)
            % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
        else
            % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
        end
        [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
        traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
        mult_factor = ones(size(traceDir));
        mult_factor(traceDir<0) = 1;
        mult_factor(traceDir>=0) = -1;
        traceND = traceDir + 90*mult_factor;    % convert traceDir to traceND
        traceSF = abs_schmid_factor(nss+1:nss+ntwin,2);
        
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        x_local = X(indR_min:indR_max, indC_min:indC_max);
        y_local = Y(indR_min:indR_max, indC_min:indC_max);
        
        % initialize for each grain (iS)
        twinMapLocal = zeros(size(ID_local));
        sfMapLocal = zeros(size(ID_local));
        
        twinMapCell_cluster = [];
        sfMapCell_cluster = [];
        
        % Used to be a function from this point on ...
        clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
        clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
        clusterNumMapC(clusterNumMapC~=iC) = 0;
        % If cluster not too close to gb
        if struCell{iE}(iS).cToGbDist(iC,end) < 50*0
            clusterNumMapC(clusterNumMapC==iC) = 0;
        end
        % If thinTF==1
        clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
        
        % Modify the 'cTrueTwin' field, modification is based on manual input    
        struCell{iE}(iS).cTrueTwin(iC,:) = activeSS;
        
        % If there are more than one active slip/twin systems, we should seperate them:
        switch sum(activeSS)
            case 0
                fragments = zeros(size(clusterNumMapC));
                
            case 1     % If try always run fit r2 procedure, can disable this part by using case '-1'
                ind = find(activeSS);
                fragments = (nss+ind) * ones(size(clusterNumMapC));
                fragments(clusterNumMapC==0) = 0;
                
            otherwise
                % (11) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
                % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
                
                % (11.1) break skeleton into small branches
                skl = clusterNumMapT;
                branchPoints = bwmorph(skl, 'branchpoints');
                
                % If heavyClean = 1;
                branchPoints = imdilate(branchPoints, ones(3));
                branch = skl - branchPoints;
                
                % (11.2) assign an ID to each skeleton branch
                branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
                branchNumbered(~branch) = 0;
                % and get the unique ID of the branches
                uniqueBranchNum = unique(branchNumbered(:));
                uniqueBranchNum(uniqueBranchNum==0)=[];
                % [illustrate] skeleton branches
                % myplotm(mod(branchNumbered,5)+logical(branchNumbered));
                
                % (11.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.
                % Assign the ts/ss ID to the branches, which can be considered as grouped.
                branchGrouped = zeros(size(branchNumbered));
                for ib = 1:length(uniqueBranchNum)
                    model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                    
                    branchND = atand(-1/model.Coefficients.Estimate(2));
                    dAngle = abs(traceND - branchND);
                    dAngle(~activeSS) = inf;
                    [~,ind] = min(dAngle);
                    branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;
                end
                
                % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
                [~,fragments] = city_block(branchGrouped);
                fragments(clusterNumMapC==0) = 0;
                
                % [illustrate] the fragments
                % myplot(fragments, branch); caxis([18,24]);
        end
        
        
        % Only update iE level.  Other levels were updated iteratively.
        if ~isempty(fragments)
            twinMapCell_cluster = fragments;
            
            sfMap = zeros(size(fragments));
            for it = 1:ntwin
                sfMap(fragments==it+nss) = traceSF(it);
            end
            sfMapCell_cluster = sfMap;
        end
        
        
        
        if ~isempty(twinMapCell_cluster)
            twinMapLocal = twinMapLocal + twinMapCell_cluster;
            sfMapLocal = sfMapLocal + sfMapCell_cluster;
        end
        
        
        % update. First clean old map, then add new map.
        map_local = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);  % (1) Cut a squared map from big map
        trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;          % (2) Eliminate this sqaured region from the big map
        map_local(clusterNumMapC==iC) = 0;                              % (3) clean the cluster (instead of grain) area on the cut map
        map_local = map_local + twinMapLocal;                           % (4) update the grain area on the cut map
        trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;  % (5) Add the modified cut map to big map
        
        map_local = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;
        map_local(clusterNumMapC==iC) = 0;
        map_local = map_local + sfMapLocal;
        sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;
        
    end
end

warning('on','all');


% Count tVol again
for iE = iE_start:iE_stop
    % use the updated trueTwinMapCell
    
    [idic,IA,IC] = unique(ID(:)*1000+trueTwinMapCell{iE}(:));
    idicCount = histcounts(ID(:)*1000+trueTwinMapCell{iE}(:),[idic(:);idic(end)+1]);
    
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        
        for iTwin = 1:length(struCell{iE}(iS).tSF)
            if struCell{iE}(iS).tSF(iTwin) < SF_th
                struCell{iE}(iS).cTrueTwin(:,iTwin) = 0;
            end
            ind = find(idic==ID_current*1000+struCell{iE}(iS).tLabel(iTwin));
            if ~isempty(ind)
                struCell{iE}(iS).tVol(iTwin) = idicCount(ind);
            end
        end
    end
    
end

%% show the relabeled, updated map
close all;
for iE = 2:5
    myplot(trueTwinMapCell{iE},boundaryTFB); caxis([18 24]);
end

%% save again
timeStr = datestr(now,'yyyymmdd_HHMM');
save([timeStr,'_relabeled_result.mat'],'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote','-v7.3');







