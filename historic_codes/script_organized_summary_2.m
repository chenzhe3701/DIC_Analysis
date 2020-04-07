% script for analyzing twin-grain boundary intersection
%
% 2019-06-13 note
% Check [strain] distribution at [dist_to_gb]
% criterion: involved_twinned, involved_not_twinned, not_involved
% This might NOT be a good way to summarize.  
% Also, this only look at the current strain level

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
% Make unique grain boundary map, and a list of the unique grain boundaries
[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:)); 
uniqueBoundaryList(uniqueBoundaryList==0) = [];
clear boundaryID neighborID;

%% (2) Check [strain] distribution at [dist_to_gb]
% myplot(distMap);  % it seems enough to just consider distance with 150 data points to grian boundary

for iE_select = iE_start:iE_stop
    
    ind = (tb_iE <= iE_select);
    activity_this_iE = unique([tb_gbNum(ind), tb_gNum(ind), tb_tsNum(ind)],'rows');
    
    bN_involved_twinned = [];
    gN_involved_twinned = [];
    bN_involved_not_twinned = [];
    gN_involved_not_twinned = [];
    bN_not_involved = [];
    gN_not_involved = [];
    
    for ii = 1:length(uniqueBoundaryList)
        this_gb = uniqueBoundaryList(ii);
        
        if ~ismember(this_gb, activity_this_iE(:,1))    % (3) This gb is not a twin gb at this iE
            gs = [floor(this_gb/10000), mod(this_gb,10000)];
            g1 = gs(1);
            g2 = gs(2);
            bN_not_involved = [bN_not_involved; this_gb; this_gb];
            gN_not_involved = [gN_not_involved; g1; g2];
        else            
            ind = activity_this_iE(:,1)==this_gb;
            activity_this_iE_gb = activity_this_iE(ind,:);
            for jj = 1:size(activity_this_iE_gb,1)
                activity = activity_this_iE_gb(jj,:);
                gs = [floor(activity(1)/10000), mod(activity(1),10000)];
                gb = activity(1);
                g_twinned = activity(2);
                g_neighbor = gs(~ismember(gs,g_twinned));
                
                % (1) Record those involved_twinned
                bN_involved_twinned = [bN_involved_twinned; gb];
                gN_involved_twinned = [gN_involved_twinned; g_twinned];
                % (2) If neighbor grain is not twinned, record it as involved_not_twinned   
                if ~ismember([gb, g_neighbor], activity_this_iE_gb(:,1:2),'rows')
                    bN_involved_not_twinned = [bN_involved_not_twinned; gb];
                    gN_involved_not_twinned = [gN_involved_not_twinned; g_neighbor];
                end

            end
        end
        
    end
    
    % (AAA) This strain level
%     ind = (tb_iE <= iE_select);   % Note: <= for accumulative ---------------------------------------------------------------------------
%     tbN = unique(tb_gbNum(ind));   % twin_gb_number at this iE_select
%     nonTwinnedGbN = uniqueBoundaryList(~ismember(uniqueBoundaryList,tbN));
%     gN_involved = unique([mod(tbN,10000); floor(tbN/10000)]);  % all grain involved
%     gN_involved_twinned = unique(tb_gNum(ind));     % twinned_grain_number at this iE_select
%     gN_involved_not_twinned = gN_involved(~ismember(gN_involved,gN_involved_twinned));  % all grains involved in twinned, the grain may just be on one side of twinGB, but no have a twin
%     gN_not_involved = gID(~ismember(gID,gN_involved));
    
    % (1) distance from twin boundary, grain is invovled_twinned
    distMap_involved_twinned = distance_from_boundary_in_grain(ID, bN_involved_twinned, gN_involved_twinned);
    
    % (2) distance from twin boundary, grain is invovled_not_twinned
    distMap_involved_not_twinned = distance_from_boundary_in_grain(ID, bN_involved_not_twinned, gN_involved_not_twinned);
    
    % (3) distance from non-twinned boundary, grain is not_involved
    distMap_not_involved = distance_from_boundary_in_grain(ID, bN_not_involved, gN_not_involved);

    % (4) emap
    eMap = calculate_effective_strain(strainFile{iE_select}.exx, strainFile{iE_select}.exy, strainFile{iE_select}.eyy);
    
    % This is from the maps, for the quantiles
    q_involved_not_twinned = [];
    q_involved_twinned = [];
    q_not_involved = [];
    
    NN = max([distMap_involved_twinned(:);distMap_involved_not_twinned(:);distMap_not_involved(:)])
    NN = 250;
    x_dist = (1:NN)';
    qs_target = [0.0014, 0.0227, 0.1587, 0.5, 0.8414, 0.9772, 0.9987];
   
    q_involved_twinned = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_involved_twinned(:)==x), qs_target), nanmean(eMap(distMap_involved_twinned(:)==x))], x_dist, 'uniformoutput', 0)))';
    q_involved_not_twinned = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_involved_not_twinned(:)==x), qs_target), nanmean(eMap(distMap_involved_not_twinned(:)==x))], x_dist, 'uniformoutput', 0)))';
    q_not_involved = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_not_involved(:)==x), qs_target), nanmean(eMap(distMap_not_involved(:)==x))], x_dist, 'uniformoutput', 0)))';

    save(['2_quants_wrt_gb_',num2str(iE_select),'.mat'], 'x_dist','q_involved_not_twinned','q_involved_twinned','q_not_involved',...
        'distMap_involved_twinned','distMap_involved_not_twinned','distMap_not_involved');
end

%% (2.1) After getting data, select an iE_select, plot the distance masks/maps  
% close all;
iE_select = 5;
load(['2_quants_wrt_gb_',num2str(iE_select),'.mat']);
um_per_dp = 5*360/4096;    % micron per data point, ~0.43
% x_dist = x_dist * um_per_dp;
myplot(distMap_involved_twinned,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_involved_not_twinned,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_not_involved,boundaryTFB + logical(trueTwinMapCell{iE_select}));

%% for the [three different regions], plot quantiles, or mean values, etc
ith = 7;
figure; hold on;
ndp_to_plot = length(x_dist);    % or 150
% ndp_to_plot = 150;
plot(x_dist(1:ndp_to_plot),q_involved_twinned(ith,1:ndp_to_plot),'-r');
plot(x_dist(1:ndp_to_plot),q_involved_not_twinned(ith,1:ndp_to_plot),'-b');
plot(x_dist(1:ndp_to_plot),q_not_involved(ith,1:ndp_to_plot),'-k');
legend({'Involved twinned','Involved not twinned','Not involved'});
set(gca,'fontsize',16);
xlabel('Distance to grain boundary');
ylabel('Effective strain');
title(['iE=',num2str(iE_select)],'fontweight','normal');


%% for a [selected region, e.g., involved_twinned], plot median and percentile
% close all;
figure; hold on;
ndp_to_plot = find(isnan(sum(q_involved_twinned,1)),1,'first')-1
ndp_to_plot = 200;

plot(x_dist(1:ndp_to_plot),q_involved_twinned(3,1:ndp_to_plot),'-r','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_involved_twinned(4,1:ndp_to_plot),'-k','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_involved_twinned(5,1:ndp_to_plot),'-b','LineWidth',2);
legend('16 percentile', 'Median', '84 percentile');
xlabel('Distance to boundary');
set(gca,'fontsize',16);
title(['iE = ',num2str(iE_select)]);
% xx = [x_dist(1:ndp_to_plot);flipud(x_dist(1:ndp_to_plot))]';
% yy = [q_involved_twinned(3,1:ndp_to_plot), fliplr(q_involved_twinned(5,1:ndp_to_plot))];
% fill(xx,yy, 'r','FaceAlpha',0.4);







