% script for analyzing twin-grain boundary intersection
%
% chenzhe, 2019-06-24 note
% summarize the strain distribution (mean, quantile) within [distance] to [gb] at [type of grains] e.g.,  
% [gb, grain - nontwin side]  twin-induced-twin activity, old twin side  
% [gb, grain - twin side] twin-induced-twin activity, new twin side
% [gb, grain - coactivated] coactivated activity
% [gb, grain -slip side] slip-induced activity, slip side 
% [gb, grain -twin side] slip-induced activity, twin side
% [gb, grain] where no twin activity is invovlde   
%
% For each grain, get one data, then use histogram to summarize.

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
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\','select the results where twin identification was based on trace dir and strain');

[twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\p\m\DIC_Analysis\20190209_1933_twin_at_boundary_result.mat','select the results for twin-grain boundary intersection');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
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

%% (0) load data, using SF threshold values to assign active twin system, and make maps
% Load cluster number maps (cleaned).
clusterNumberMapCell = cell(1,length(STOP)-1);
for iE = []%iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
% load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');
load(fullfile(confirmedLabelPath,confirmedLabelFile),'trueTwinMapCell');
%% load previous twin_gb interaction result
load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile));

%% (5) select grains that will twin, look at strain in previous strain level (grain-wise summary)
for iE_select = 4%iE_start+1:iE_stop
    
    load(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat'], 'distMap_twin_twin_n_side','distMap_twin_twin_t_side','distMap_coactivation','distMap_slip_twin_n_side','distMap_slip_twin_t_side','distMap_not_involved');
    
    % () Find activities [gbN, gN, tsN] at adjacent [iEs]
    ind = (tb_iE <= iE_select-1);
    activity_previous_iE = unique([tb_gbNum(ind), tb_gNum(ind), tb_tsNum(ind)],'rows');
    
    ind = (tb_iE <= iE_select);
    activity_this_iE = unique([tb_gbNum(ind), tb_gNum(ind), tb_tsNum(ind)],'rows');
    
    % (5) emap
    eMap = calculate_effective_strain(strainFile{iE_select-1}.exx, strainFile{iE_select-1}.exy, strainFile{iE_select-1}.eyy);
    
    % pn contains data pair [grain boundary number, grain number of grain that satisfies some criterion: (newly activated events at this strain)   
    p1 = [];    % [gb, grain - nontwin side]  twin-induced-twin activity, old twin side  
    p2 = [];    % [gb, grain - twin side] twin-induced-twin activity, new twin side
    p3 = [];    % [gb, grain - coactivated] coactivated activity
    p4 = [];    % [gb, grain -slip side] slip-induced activity, slip side 
    p5 = [];    % [gb, grain -twin side] slip-induced activity, twin side
    p6 = [];    % [gb, grain] where no twin activity is invovlde   
    for ii = 1:size(activity_this_iE,1)
        activity = activity_this_iE(ii,:);
        % (1) If activity is new? --> analyze newly occured activity
        if ~ismember(activity,activity_previous_iE,'rows')
            gs = [floor(activity(1)/10000), mod(activity(1),10000)];
            gb = activity(1);
            g_twinned = activity(2);
            g_neighbor = gs(~ismember(gs,g_twinned));
                        
            if ismember([gb,g_neighbor], activity_previous_iE(:,1:2))   % (1) is neighbor twinned at previous iE? --> twin induced twin activity.
                p1 = [p1; gb, g_neighbor];
                p2 = [p2; gb, g_twinned];
            elseif ismember([gb,g_neighbor], activity_this_iE(:,1:2))   % (2) is neighbor twinned (so just twinned) at this iE? --> coactivated twin activity.
                p3 = [p3; gb, g_twinned; gb, g_neighbor];
            else    % (3) neighbor not twinned yet, slip induced twin activity
                p4 = [p4; gb, g_neighbor];
                p5 = [p5; gb, g_twinned];
            end
        end
    end
    for ii = 1:length(uniqueBoundaryList)
        gb = uniqueBoundaryList(ii);
        if ~ismember(gb, activity_this_iE(:,1))    % (3) This gb is not a twin gb at this iE
            gs = [floor(gb/10000), mod(gb,10000)];
            g1 = gs(1);
            g2 = gs(2);
            p6 = [p6; gb, g1; gb, g2];
        end
    end

    p1 = unique(p1,'rows');
    p2 = unique(p2,'rows');
    p3 = unique(p3,'rows');
    p4 = unique(p4,'rows');
    p5 = unique(p5,'rows');
    p6 = unique(p6,'rows');
  
    th = 100;
    d4 = [];    % 'd' for data.  Summarize mean/quantile of strain 
    d44 = [];   
    for ii = 1:size(p4,1)
        ii
        ID_interest = p4(ii,2);
        ind_local = ismember(ID, ID_interest); %ismember(ID, [ID_current,ID_neighbor]);
        
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        distMap_local = distMap_slip_twin_n_side(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);
        
        data = eMap_local(distMap_local<th);
        qs_target = [0.0014, 0.0227, 0.1587, 0.5, 0.8414, 0.9772, 0.9987];
        
        d4 = [d4; nanmean(data)];
        d44 = [d44, quantile(data,0.5)];
    end
    
    d5 = [];
    d55 = [];
    for ii = 1:size(p5,1)
        ii
        ID_interest = p5(ii,2);
        ind_local = ismember(ID, ID_interest); %ismember(ID, [ID_current,ID_neighbor]);
        
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        distMap_local = distMap_slip_twin_t_side(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);
        
        data = eMap_local(distMap_local<th);
        qs_target = [0.0014, 0.0227, 0.1587, 0.5, 0.8414, 0.9772, 0.9987];
        
        d5 = [d5; nanmean(data)];
        d55 = [d55, quantile(data,0.5)];
    end

    figure; histogram(d44,0:0.001:0.06);  set(gca,'xlim',[0 0.06])
    figure; histogram(d55,0:0.001:0.06);  set(gca,'xlim',[0 0.06])
end





