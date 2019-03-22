% script for analyzing twin-grain boundary intersection

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

%% (3) select grains that will twin, look at strain in previous strain level
for iE_select = iE_start+1:iE_stop
    
    % () Find activities [gbN, gN, tsN] at adjacent [iEs]
    ind = (tb_iE <= iE_select-1);
    activity_previous_iE = unique([tb_gbNum(ind), tb_gNum(ind), tb_tsNum(ind)],'rows');
    
    ind = (tb_iE <= iE_select);
    activity_this_iE = unique([tb_gbNum(ind), tb_gNum(ind), tb_tsNum(ind)],'rows');
    
    bN_twin_twin = [];
    gN_twin_twin_n_side = [];
    gN_twin_twin_t_side = [];
    bN_coactivation = [];
    gN_coactivation = [];
    bN_slip_twin = [];
    gN_slip_twin_n_side = [];
    gN_slip_twin_t_side = [];
    bN_not_involved = [];
    gN_not_involved = [];
    
    for ii = 1:size(activity_this_iE,1)
        activity = activity_this_iE(ii,:);
        % (1) If activity is new?
        if ~ismember(activity,activity_previous_iE,'rows')
            gs = [floor(activity(1)/10000), mod(activity(1),10000)];
            gb = activity(1);
            g_twinned = activity(2);
            g_neighbor = gs(~ismember(gs,g_twinned));
                        
            if ismember([gb,g_neighbor], activity_previous_iE(:,1:2))   % (1) is neighbor twinned at previous iE? --> twin induced twin activity.
                bN_twin_twin = [bN_twin_twin; gb];
                gN_twin_twin_n_side = [gN_twin_twin_n_side; g_neighbor];                
                gN_twin_twin_t_side = [gN_twin_twin_t_side; g_twinned];
            elseif ismember([gb,g_neighbor], activity_this_iE(:,1:2))   % (2) is neighbor twinned (so just twinned) at this iE? --> coactivated twin activity.
                bN_coactivation = [bN_coactivation; gb; gb];
                gN_coactivation = [gN_coactivation; g_neighbor; g_twinned]; % both are involved in coactivation
            else    % (3) neighbor not twinned yet, slip induced twin activity
                bN_slip_twin = [bN_slip_twin; gb];
                gN_slip_twin_n_side = [gN_slip_twin_n_side; g_neighbor];
                gN_slip_twin_t_side = [gN_slip_twin_t_side; g_twinned];
            end
        end
    end
    for ii = 1:length(uniqueBoundaryList)
        this_gb = uniqueBoundaryList(ii);
        if ~ismember(this_gb, activity_this_iE(:,1))    % (3) This gb is not a twin gb at this iE
            gs = [floor(this_gb/10000), mod(this_gb,10000)];
            g1 = gs(1);
            g2 = gs(2);
            bN_not_involved = [bN_not_involved; this_gb; this_gb];
            gN_not_involved = [gN_not_involved; g1; g2];
        end
    end

    % (1) distance from twin boundary, grain is invovled_twinned
    distMap_twin_twin_n_side = distance_from_boundary_in_grain(ID, bN_twin_twin, gN_twin_twin_n_side);
    distMap_twin_twin_t_side = distance_from_boundary_in_grain(ID, bN_twin_twin, gN_twin_twin_t_side);
    % (2) distance from twin boundary, grain is invovled_not_twinned
    distMap_coactivation = distance_from_boundary_in_grain(ID, bN_coactivation, gN_coactivation);
    
    % (3) distance from non-twinned boundary, grain is not_involved
    distMap_slip_twin_n_side = distance_from_boundary_in_grain(ID, bN_slip_twin, gN_slip_twin_n_side);
    distMap_slip_twin_t_side = distance_from_boundary_in_grain(ID, bN_slip_twin, gN_slip_twin_t_side);
    
    % (4) distance from non-twinned boundary, grain is not_involved
    distMap_not_involved = distance_from_boundary_in_grain(ID, bN_not_involved, gN_not_involved);
    
    % (5) emap
    eMap = calculate_effective_strain(strainFile{iE_select-1}.exx, strainFile{iE_select-1}.exy, strainFile{iE_select-1}.eyy);
    
    % This is from the maps, for the quantiles
    NN = max([distMap_twin_twin_n_side(:);distMap_twin_twin_t_side(:);distMap_coactivation(:);distMap_slip_twin_n_side(:);distMap_slip_twin_t_side(:)])
    NN = 250;
    x_dist = (1:NN)';
    qs_target = [0.0014, 0.0227, 0.1587, 0.5, 0.8414, 0.9772, 0.9987];
   
    q_twin_twin_n_side = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_twin_twin_n_side(:)==x), qs_target), nanmean(eMap(distMap_twin_twin_n_side(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 1');
    q_twin_twin_t_side = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_twin_twin_t_side(:)==x), qs_target), nanmean(eMap(distMap_twin_twin_t_side(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 2');
    q_coactivation = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_coactivation(:)==x), qs_target), nanmean(eMap(distMap_coactivation(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 3');
    q_slip_twin_n_side = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_twin_n_side(:)==x), qs_target), nanmean(eMap(distMap_slip_twin_n_side(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 4');
    q_slip_twin_t_side = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_twin_t_side(:)==x), qs_target), nanmean(eMap(distMap_slip_twin_t_side(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 5');
    q_not_involved = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_not_involved(:)==x), qs_target), nanmean(eMap(distMap_not_involved(:)==x))], x_dist, 'uniformoutput', 0)))';
    disp('sets of quants calculated: 6');
    
    save(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat'], 'x_dist','q_twin_twin_n_side','q_twin_twin_t_side','q_coactivation','q_slip_twin_n_side','q_slip_twin_t_side','q_not_involved',...
        'distMap_twin_twin_n_side','distMap_twin_twin_t_side','distMap_coactivation','distMap_slip_twin_n_side','distMap_slip_twin_t_side','distMap_not_involved');
%     save(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat'], 'q_not_involved','distMap_not_involved','-append');
end

%% (3.1) After getting data, select an iE_select, plot the distance masks/maps  
% close all;
iE_select = 3;
load(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat']);
um_per_dp = 5*360/4096;    % micron per data point, ~0.43
myplot(logical(trueTwinMapCell{iE_select-1}), boundaryTFB);
myplot(strainFile{iE_select-1}.exx, boundaryTFB);
myplot(logical(trueTwinMapCell{iE_select}), boundaryTFB);
myplot(strainFile{iE_select}.exx, boundaryTFB);

myplot(distMap_twin_twin_n_side,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_twin_twin_t_side,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_coactivation,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_slip_twin_n_side,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_slip_twin_t_side,boundaryTFB + logical(trueTwinMapCell{iE_select}));
myplot(distMap_not_involved,boundaryTFB + logical(trueTwinMapCell{iE_select}));
% caxis([0 200]); set(gca,'xlim',[1600,2600],'ylim',[400 1000]); % to crop for summary   
%% for the [three different regions], plot quantiles, or mean values, etc
iE_select = 5;
load(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat']);
ith = 4;
figure; hold on;
ndp_to_plot = length(x_dist);    % or 150
% ndp_to_plot = 150;
plot(x_dist(1:ndp_to_plot),q_twin_twin_n_side(ith,1:ndp_to_plot),'-r','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_twin_twin_t_side(ith,1:ndp_to_plot),'--r','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_coactivation(ith,1:ndp_to_plot),'-k','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_slip_twin_n_side(ith,1:ndp_to_plot),'-b','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_slip_twin_t_side(ith,1:ndp_to_plot),'--b','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q_not_involved(ith,1:ndp_to_plot),'-','color',[0 0.6 0],'LineWidth',2);
legend({'twin-twin old-twin-side','twin-twin new-twin-side','co-found','slip-twin slip-side','slip-twin new-twin-side','non-involved'});
set(gca,'fontsize',16);
xlabel('Distance to grain boundary');
ylabel('Effective strain');
title(['iE=',num2str(iE_select)],'fontweight','normal');
set(gca,'ylim',[0.015 0.04]);

%% for a [selected region, e.g., involved_twinned], plot median and percentile
close all;
iE_select = 5;
load(['tracked_quants_wrt_gb_',num2str(iE_select),'.mat']);

% ndp_to_plot = find(isnan(sum(q_coactivation,1)),1,'first')-1
ndp_to_plot = length(x_dist); 
ndp_to_plot = 200;

for toPlot = 1:5
switch toPlot
    case 1
        q = q_twin_twin_n_side;
        title_string = ['twin-twin non-twin-side, iE = ',num2str(iE_select)];
    case 2
        q = q_twin_twin_t_side;
        title_string = ['twin-twin twin-side, iE = ',num2str(iE_select)];
    case 3
        q = q_coactivation;
        title_string = ['co-found, iE = ',num2str(iE_select)];
    case 4
        q = q_slip_twin_n_side;
        title_string = ['slip-twin slip-side, iE = ',num2str(iE_select)];
    case 5
        q = q_slip_twin_t_side;
        title_string = ['slip-twin twin-side, iE = ',num2str(iE_select)];
    case 6
        q = q_not_involved;
        title_string = ['slip-twin twin-side, iE = ',num2str(iE_select)];
    otherwise   
end

legend_string = {'2 percentile','16 percentile', 'Median', '84 percentile','98 percentile','mean'}; % [0.1, 2, 16, 50, 84, 98, 99.9] percentile 

figure; hold on;
plot(x_dist(1:ndp_to_plot),q(2,1:ndp_to_plot),'--r','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q(3,1:ndp_to_plot),'-r','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q(4,1:ndp_to_plot),'-k','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q(5,1:ndp_to_plot),'-b','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q(6,1:ndp_to_plot),'--b','LineWidth',2);
plot(x_dist(1:ndp_to_plot),q(8,1:ndp_to_plot),'-g','LineWidth',2);
legend(legend_string);  
% legend('16 percentile', 'Median', '84 percentile','mean');
xlabel('Distance to boundary (# data points)');
set(gca,'fontsize',16);
title(title_string,'fontweight','normal');
% xx = [x_dist(1:ndp_to_plot);flipud(x_dist(1:ndp_to_plot))]';
% yy = [q_involved_twinned(3,1:ndp_to_plot), fliplr(q_involved_twinned(5,1:ndp_to_plot))];
% fill(xx,yy, 'r','FaceAlpha',0.4);
set(gca,'ylim',[0 0.1]);

end




