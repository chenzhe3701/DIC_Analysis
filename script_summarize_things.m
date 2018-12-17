
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
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\','select the results where twin identification was based on trace dir and strain');

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
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');

%% can show map
close all;
for iE = 2:5
    myplot(trueTwinMapCell{iE},boundaryTFB); caxis([18 24]);
end


%% Summarize, e.g., grain with twins, without twins, and their distribution on IPF/pole figure
%% (s0) The SF distribution of all grains w.r.t grain size. (i.e., does larger grain tend to have larger/smaller SF? should not, but was asked to check ...)
% (1) using box plot to check
% (2) plot on IPF map to check
close all;
um_per_dp = 5*360/4096;    % micron per data point, ~0.43
stru = struCell{2};
eulers = [];
sf = [];
gd = [];
gs = [];

for iS = 1:length(stru)
    ID_current = stru(iS).gID;
    ind = find(gID == ID_current);
    eulers(iS,:) = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
    sf(iS) = max(stru(iS).tSF);
    gd(iS) = sqrt(struCell{iE}(iS).gVol/pi*4) * um_per_dp;
    gs(iS) = struCell{iE}(iS).gVol*um_per_dp.^2;
end

edges_d = [0:20:20*15,400];
edges_a = [0:3000:3000*15, 125000];

figure; histogram(gd,edges_d);ylabel('counts');
label_d = discretize(gd,edges_d);
hold on;
yyaxis right; boxplot(sf, label_d, 'positions', 10:20:20*16); 
for ii = 1:length(edges_d)-1
    labels_d{ii} = [num2str(edges_d(ii)),'-',num2str(edges_d(ii+1))];
end
set(gca,'xticklabels',labels_d,'xticklabelrotation',45);
ylabel('Schmid Factor (m)');
xlabel('Grain diameter (um)');

for ii = unique(label_d)
    plot_on_IPF(eulers(label_d==ii,:),[0 0 0],[0 0 0],[1 0 0],19:24,[-1 0 0; 0 0 0; 0 0 0],'Mg','Twin');
    title([labels_d{ii},'um']);
end

figure; histogram(gs,edges_a);ylabel('counts');
label_a = discretize(gs,edges_a);
hold on;
yyaxis right; boxplot(sf, label_a, 'positions', 1500:3000:3000*16);
for ii = 1:length(edges_d)-1
    labels_a{ii} = [num2str(edges_a(ii)),'-',num2str(edges_a(ii+1))];
end
set(gca,'xticklabels',labels_a,'xticklabelrotation',45);
ylabel('Schmid Factor (m)');
xlabel('Grain size (um^2)');

%% (s1) as discussed on Nov-14, check [grains_twinned, grains-nontwinned] as a function of [grain_size], -- large vs small grains.

close all;
for iE = iE_start:iE_stop
    um_per_dp = 5*360/4096;    % micron per data point, ~0.43
    gs_t = [];
    gs_nt = [];
    
    for iS = 1:length(struCell{iE})
        if any(sum(struCell{iE}(iS).cTrueTwin,1))
            gs_t = [gs_t, struCell{iE}(iS).gVol*um_per_dp.^2];
        else
            gs_nt = [gs_nt, struCell{iE}(iS).gVol*um_per_dp.^2];
        end
    end
    
    edges = 0:1000:20000;
    figure;
    histogram([gs_nt,gs_t],edges);
    hold on;
    histogram(gs_t,edges);
    legend('Total', 'Twinned','location','west');
    xlabel('grain size');
    ylabel('Counts');
    set(gca,'ylim',[0 120]);
    title(['iE = ',num2str(iE)]);
    yyaxis right;
    set(gca,'ycolor','r')
    plot(edges(1:end-1)+(edges(2)-edges(1))/2,histcounts(gs_t,edges)./(histcounts(gs_t,edges)+histcounts(gs_nt,edges)),'-ro');
    ylabel('pct grain twinned');
    set(gca,'ylim',[0 0.4]);
end

%%
close all;
for iE = iE_start:iE_stop
    um_per_dp = 5*360/4096;    % micron per data point, ~0.43
    gs_t = [];
    sf_t = [];
    gs_nt = [];
    sf_nt = [];
    
    ic_t = 1;
    ic_nt = 1;
    for iS = 1:length(struCell{iE})
        if any(sum(struCell{iE}(iS).cTrueTwin,1))
            gs_t(ic_t) = struCell{iE}(iS).gVol*um_per_dp.^2;
            sf_t(ic_t) = mean(struCell{iE}(iS).tSF(sum(struCell{iE}(iS).cTrueTwin,1)>0));
            ic_t = ic_t + 1;
        else
            gs_nt(ic_nt) = struCell{iE}(iS).gVol*um_per_dp.^2;
            sf_nt(ic_nt) = max(struCell{iE}(iS).tSF);
            ic_nt = ic_nt + 1;
        end
    end
    
%     edges_d = [0:20:20*15,400];
    edges_a = [0:3000:3000*15, 125000];

    bn_t = discretize(gs_t,edges_a);
    bn_nt = discretize(gs_nt,edges_a);
    
    figure;
    histogram([gs_nt,gs_t],edges_a);
    hold on;
    histogram(gs_t,edges_a);

    xlabel('grain size (um^2)');
    ylabel('Counts');
    title(['iE = ',num2str(iE)]);
    yyaxis right;
    set(gca,'ycolor','r')
    plot(edges_a(1:end-1)+(edges_a(2)-edges_a(1))/2,histcounts(gs_t,edges_a)./(histcounts(gs_t,edges_a)+histcounts(gs_nt,edges_a)),'-ro');
    ylabel('pct grain twinned');
    set(gca,'ylim',[0 0.4], 'xlim',[0,48000]);
    legend('Total', 'Twinned','Pct Grains Twinned','location','west');
    
    for ii = 1:length(edges_d)-1
        labels_a{ii} = [num2str(edges_a(ii)),'-',num2str(edges_a(ii+1))];
    end
    figure;
    boxplot(sf_t,bn_t);
    set(gca,'xticklabels',labels_a,'xticklabelrotation',45);
    title('SF Distribution in Twinned Grains');
    xlabel('grain size (um^2)'); ylabel('Schmid Factor (m)');
    
    figure;
    boxplot(sf_nt,bn_nt);
    set(gca,'xticklabels',labels_a,'xticklabelrotation',45);
    title('SF Distribution in Non-Twinned Grains');
    xlabel('grain size (um^2)'); ylabel('Schmid Factor (m)');
end





%% (s2) Only consider twinned grains, [pct_of_twinned_area] as a function of [grain_size] 
close all;
for iE = iE_start:iE_stop
    um_per_dp = 5*360/4096;    % micron per data point, ~0.43
    tVol = [];
    gVol = [];
    gSize = [];
    
    for iS = 1:length(struCell{iE})
        if any(sum(struCell{iE}(iS).cTrueTwin,1))
            tVol = [tVol, sum(struCell{iE}(iS).tVol)];
            gVol = [gVol, struCell{iE}(iS).gVol];
            if tVol(end)/gVol(end)>1
                disp(['found ID:',num2str(iS)]);
            end
            gSize = [gSize, struCell{iE}(iS).gVol * um_per_dp.^2];
        end
    end
    
    figure;
    plot(gSize,tVol./gVol,'or');
    set(gca,'xlim',[0 120000]);
    xlabel('grain size (um^2) of twinned grain');
    ylabel('twin area pct');
    title(['iE = ',num2str(iE)]);
end


%% (s3) twinned grain]


%% (2) For a twinned grain, what is the SF of that twin.  For a non-twinned grain, what is the maximum SF.
close all;

% investigate the grain size, get an estimate, and eliminate the small grains  
gSizeMap = zeros(size(ID));
gSizes = zeros(1,length(struCell{5}));
for ii = 1:length(struCell{5})
   gSizes(ii) = struCell{5}(ii).gVol;
   gSizeMap(ID==struCell{5}(ii).gID) = struCell{5}(ii).gVol;
end
figure; histogram(gSizes);
myplotc_low(gSizeMap);  % choose 3000?

% using this gSizeMap, we can quantile(gSizeMap(:),pct) to get the gSize below which accounts for pct of the total area 
quantile(gSizeMap(:),[0.05, 0.4, 0.5, 0.6, 0.95])

save('gSizeMap.mat','gSizeMap','gSizes');
%% (2.1) summarize. from previous section, for WE43_T6, there is not a minimum grain size to consider. But it is possible for other sample  
for iE = iE_start:iE_stop
    mTwinned = [];
    mNonTwinned = [];

    for iS = 1:length(struCell{iE})
        t = sum(struCell{iE}(iS).cTrueTwin,1);
        ind = find(t(:));
        if isempty(ind)
            mNonTwinned = [mNonTwinned; max(struCell{iE}(iS).tSF)];
        else
            t_SF = struCell{iE}(iS).tSF(ind);
            mTwinned = [mTwinned; t_SF(:)];
        end
    end
    
    figure;
    histogram(mTwinned,-0.5:0.05:0.5);
    hold on;
    histogram(mNonTwinned,-0.5:0.05:0.5);
    legend('Twinned','Not Twinned','location','west');
    xlabel('m');
    ylabel('Counts');
    set(gca,'ylim',[0 350]);
    title(['iE = ',num2str(iE)]);    
    
end


%% (3) for twinned grains, what is the rank of the twin system ?

%% (4) exx_of_twin_area vs. exx_of_nontwinned_area.  And total u accommodated by twin vs nontwinned area.
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
for iE = iE_start:iE_stop
    
    exx = strainFile{iE}.exx;
    exxTwinned = exx(trueTwinMapCell{iE}>nss);
    exxNonTwinned = exx(trueTwinMapCell{iE}<=nss);
    exxTwinned(isnan(exxTwinned)) = [];
    exxNonTwinned(isnan(exxNonTwinned)) = [];
    
    qt1 = quantile(exxTwinned,[0.01, 0.99]);
    qt2 = quantile(exxNonTwinned,[0.01, 0.99]);
    xi = linspace(min(qt1(1),qt2(1)), max(qt1(2),qt2(2)), 101);
    epdf_twin = ksdensity(exxTwinned,xi);
    epdf_nontwin = ksdensity(exxNonTwinned,xi);
    figure; hold on;
    plot(xi,epdf_twin,'s-');
    plot(xi,epdf_nontwin,'o-');
    legend('exx Twinned','exx Not Twinned','location','west');
    xlabel('exx');
    ylabel('pct');
    title(['iE = ',num2str(iE)]);

    
    
    exy = strainFile{iE}.exy;
    exyTwinned = exy(trueTwinMapCell{iE}>nss);
    exyNonTwinned = exy(trueTwinMapCell{iE}<=nss);
    exyTwinned(isnan(exyTwinned)) = [];
    exyNonTwinned(isnan(exyNonTwinned)) = [];
    
    qt1 = quantile(exyTwinned,[0.01, 0.99]);
    qt2 = quantile(exyNonTwinned,[0.01, 0.99]);
    xi = linspace(min(qt1(1),qt2(1)), max(qt1(2),qt2(2)), 101);
    epdf_twin = ksdensity(exyTwinned,xi);
    epdf_nontwin = ksdensity(exyNonTwinned,xi);
    figure; hold on;
    plot(xi,epdf_twin,'s-');
    plot(xi,epdf_nontwin,'o-');
    legend('exy Twinned','exy Not Twinned','location','west');
    xlabel('exy');
    ylabel('pct');
    title(['iE = ',num2str(iE)]);
    
    
    
    eyy = strainFile{iE}.eyy;
    eyyTwinned = eyy(trueTwinMapCell{iE}>nss);
    eyyNonTwinned = eyy(trueTwinMapCell{iE}<=nss);
    eyyTwinned(isnan(eyyTwinned)) = [];
    eyyNonTwinned(isnan(eyyNonTwinned)) = [];
    
    qt1 = quantile(eyyTwinned,[0.01, 0.99]);
    qt2 = quantile(eyyNonTwinned,[0.01, 0.99]);
    xi = linspace(min(qt1(1),qt2(1)), max(qt1(2),qt2(2)), 101);
    epdf_twin = ksdensity(eyyTwinned,xi);
    epdf_nontwin = ksdensity(eyyNonTwinned,xi);
    figure; hold on;
    plot(xi,epdf_twin,'s-');
    plot(xi,epdf_nontwin,'o-');
    legend('eyy Twinned','eyy Not Twinned','location','west');
    xlabel('eyy');
    ylabel('pct');
    title(['iE = ',num2str(iE)]);

    
end

% calculate u, original method
for iE = iE_start:iE_stop
    exx = strainFile{iE}.exx;
    dx  = strainFile{iE}.x(1,2) - strainFile{iE}.x(1,1);
    exxTwinnedMap = exx;
    exxTwinnedMap(trueTwinMapCell{iE}==0) = 0;
    exx_along_x = nansum(exxTwinnedMap,1)./size(exxTwinnedMap,1);
    uTwinned(iE) = nansum(exx_along_x*dx);
%     uTotal = nansum(exxTwinnedMap(:)*dx);
%     uTwinned(iE) = uTotal/nSamplePointsW;
    
    exxNonTwinnedMap = exx;
    exxNonTwinnedMap(trueTwinMapCell{iE}>0) = 0;
    exx_along_x = nansum(exxNonTwinnedMap,1)./size(exxNonTwinnedMap,1);
    uNonTwinned(iE) = nansum(exx_along_x*dx);
%     uTotal = nansum(exxNonTwinnedMap(:)*dx);
%     uNonTwinned(iE) = uTotal/nSamplePointsW;
end

% to double check, approximate u range
for iE = iE_start:iE_stop
    urange(iE) = nanmean(strainFile{iE}.u(:,end-10)) - nanmean(strainFile{iE}.u(:,10));
end

figure; hold on;
plot(iE_start:iE_stop, uTwinned(iE_start:iE_stop), '-or');
plot(iE_start:iE_stop, uNonTwinned(iE_start:iE_stop), '-xb');
plot(iE_start:iE_stop, urange(iE_start:iE_stop),'-ok');
legend({'twinned','non-twinned','approximate total'});
xlabel('strain level');
ylabel('u');
title('approximate u');

figure;
plot(iE_start:iE_stop, uTwinned(iE_start:iE_stop)./uNonTwinned(iE_start:iE_stop), '-or');
legend('uTwin/uNonTwin');
xlabel('strain level');
ylabel('uTwin/uNonTwin');
title('approximate u ratio of twinned/nontwinned')




%% Is the exx strain distribution of twinned area related to grain size ? (Daly asked 2018-11-08)
for iE = iE_start:iE_stop
    clear grain_size twin_area_fraction e_avg e_std
    iCount = 1;
    exx = strainFile{iE}.exx;
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        
        % if there is active twin system
        if sum(struCell{iE}(iS).cTrueTwin(:)) > 0
            grain_size(iCount) = struCell{iE}(iS).gVol;
            
            ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
            indC_min = find(sum(ind_local, 1), 1, 'first');
            indC_max = find(sum(ind_local, 1), 1, 'last');
            indR_min = find(sum(ind_local, 2), 1, 'first');
            indR_max = find(sum(ind_local, 2), 1, 'last');
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
            
            trueTwinMap_local = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
            trueTwinMap_local(ID_local~=ID_current) = 0;
            x_local = X(indR_min:indR_max, indC_min:indC_max);
            y_local = Y(indR_min:indR_max, indC_min:indC_max);
            

            exx_local = exx(indR_min:indR_max, indC_min:indC_max);
            ind_twin = (trueTwinMap_local>0);
            
            twin_area_fraction(iCount) = sum(trueTwinMap_local(:)>0)/grain_size(iCount);
            eData = exx_local(ind_twin);
            e_avg(iCount) = nanmean(eData);
            e_std(iCount) = nanstd(eData);
            e_stats(iCount,:) = quantile(eData,[0.023, 0.159, 0.5, 0.841, 0.977]);
            
            iCount = iCount + 1;
        end
    end

%     figure;
%     scatter3(grain_size,twin_area_fraction,e_avg);
%     xlabel('grain size');
%     ylabel('twin area fraction');
%     zlabel('e avg');
%     title(['iE=',num2str(iE)]);
%     
    figure; 
    ind = grain_size > 100000;
    xd = twin_area_fraction(ind);
    yd = e_stats(ind,:);
    errorbar(xd,yd(:,3),yd(:,1),yd(:,5),'or','markeredgecolor','k','markerfacecolor','r');
    switch iE
        case {2,3,4}
            set(gca,'xlim',[0 0.8],'ylim',[-0.25,0]);
        case 5
            set(gca,'xlim',[0 0.99],'ylim',[-0.3,0]);
    end
    xlabel('twin area fraction');
    ylabel('exx distribution');
    title('grains > 100,000 data points')
    
    figure;
    xd = twin_area_fraction(~ind);
    yd = e_stats(~ind,:);
    errorbar(xd,yd(:,3),yd(:,1),yd(:,5),'ob','markeredgecolor','k','markerfacecolor','b');
    switch iE
        case 2
            set(gca,'xlim',[0 0.3],'ylim',[-0.15,0]);
        case 3
            set(gca,'xlim',[0 0.4],'ylim',[-0.25,0]);
        case 4
            set(gca,'xlim',[0 0.8],'ylim',[-0.25,0]);
        case 5
            set(gca,'xlim',[0 0.99],'ylim',[-0.3,0]);
    end
    xlabel('twin area fraction');
    ylabel('exx distribution');
    title('grains < 100,000 data points')
end
