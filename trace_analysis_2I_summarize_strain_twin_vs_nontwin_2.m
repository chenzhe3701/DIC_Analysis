
% chenzhe, 2018-5-30
% based on trace_analysis_2N_().
% reduce redundant codes, just to summarize something of interest.
% This code plot on IPF. Can use different ways. For example, 
% Size of marker = pixel number.  Color = pct of strain by twinned pixel.

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
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx','gPhi1','gPhi','gPhi2');
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

summaryPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','Choose/make a parent path for output files and figures.'),'\'];

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');

iE_of_interest = 5;

%% (0) initialize something
useStrain = 1;
makeNewMap = 0;
makeFigure = 0;
summarizeTruth = 0;

struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    struCell{iE} = stru;

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
    eeff = calculate_effective_strain(exx,exy,eyy);
    qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
    qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
    qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
    qt_eeff = quantile(eeff(:),[0.0013,0.9987]); qt_eeff(1)=min(-1,qt_eeff(1)); qt_eeff(2)=max(1,qt_eeff(2));
    ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
    exx(ind_outlier) = nan;
    exy(ind_outlier) = nan;
    eyy(ind_outlier) = nan;
    eeff(ind_outlier) = nan;
    
    EXX{iE} = exx;
    EXY{iE} = exy;
    EYY{iE} = eyy;
    EEFF{iE} = eeff;
    
end
%% select iE to analyze


% step = 0.01;
% [xMap,yMap] = meshgrid(0:step:1,0:step:1);
% twinPixelMap = zeros(size(xMap));
% nonTwinPixelMap = zeros(size(xMap));
% twinCumStrainMap = zeros(size(xMap));
% nonTwinCumStrainMap = zeros(size(xMap));

for iE = iE_start:iE_stop

    % load data for this iE
    warning('off','all');
    
    % use this as the ground truth. Note that the field 'tProb' and the map 'cnnTwinMap' might need to be processed by cnn.
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result], 'stru','clusterNumMapCleaned','trueTwinMap');
    
    n = length(stru);
    s_x = zeros(n,1);
    s_y = zeros(n,1);
    s_euler = zeros(n,3);
    s_tp = zeros(n,1);  % number of twin pixel
    s_ts = zeros(n,1);  % accumulative twin strain
    s_ap = zeros(n,1);  % number of all pixel
    s_as = zeros(n,1);  % accumulative all strain
    s_fp = zeros(n,1);  % fraction of pixels that twinned
    s_fs = zeros(n,1);  % fraction of strain accommodated by twin

    exx = EXX{iE};
    exy = EXY{iE};
    eyy = EYY{iE};
    eeff = EEFF{iE};

    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(stru)
        % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        % iS = find(gIDwithTrace == 296); % for debugging.
        ID_current = gIDwithTrace(iS);
        euler = [gPhi1(find(gID==ID_current)),gPhi(find(gID==ID_current)),gPhi2(find(gID==ID_current))];
        [posX,posY] = calculate_pos_on_IPF(euler,sampleMaterial);
%         indc = ceil((posX-xMap(1))/step + 1);
%         indr = ceil((posy-yMap(1))/step + 1);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);

        
        
        exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
        exy_local = exy(indR_min:indR_max, indC_min:indC_max);
        eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        x_local = X(indR_min:indR_max, indC_min:indC_max);
        y_local = Y(indR_min:indR_max, indC_min:indC_max);
        
        
        exx_local(ID_local~=ID_current) = nan;  % make eij anum only for the current grain
        exy_local(ID_local~=ID_current) = nan;
        eyy_local(ID_local~=ID_current) = nan;
        eeff_local = calculate_effective_strain(exx_local,exy_local,eyy_local);
        
        
        clusterNumMapLocal = clusterNumMapCleaned(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        
        % [code here] to summarize something for each cluster
        
        nCluster = length(stru(iS).cLabel);
        nTwinPixels = 0;
        nNonTwinPixels = 0;
        cumStrainTwin = 0;
        cumStrainNonTwin = 0;
        gVol = stru(iS).gVol;   % actual grain size
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            cVol = stru(iS).cVolCleaned(iCluster);
            
            if cVol>0
                % calculate the within-cluster-distance: a vector d
                indClusterLocal = (clusterNumMapLocal==cNum);
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>=0)
                    nTwinPixels = nTwinPixels + sum(indClusterLocal(:));
                    cumStrainTwin = cumStrainTwin + sum(eeff_local(indClusterLocal));
                else
                    nNonTwinPixels = nNonTwinPixels + sum(indClusterLocal(:));
                    cumStrainNonTwin = cumStrainTwin + sum(eeff_local(indClusterLocal));
                end
            end

        end
        
        s_x(iS) = posX;
        s_y(iS) = posY;
        s_euler(iS,:) = euler;
        s_tp(iS) = nTwinPixels;
        s_ts(iS) = cumStrainTwin;
        s_ap(iS) = nTwinPixels + nNonTwinPixels;
        s_as(iS) = cumStrainTwin + cumStrainNonTwin;
        
%         twinPixelMap(indr,indc) = twinPixelMap(indr,indc) + nTwinPixels;
%         nonTwinPixelMap(indr,indc) = nonTwinPixelMap(indr,indc) + nNonTwinPixels;
%         
%         twinCumStrainMap(indr,indc) = twinCumStrainMap(indr,indc) + cumStrainTwin;
%         nonTwinCumStrainMap(indr,indc) = nonTwinCumStrainMap(indr,indc) + cumStrainNonTwin;
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    
    s_fp = s_tp./s_ap;
    s_fs = s_ts./s_as;
    
%     allPixelMap = twinPixelMap + nonTwinPixelMap;
%     allCumStrainMap = twinCumStrainMap + nonTwinCumStrainMap;
%     
%     allPixelMap(allPixelMap==0) = nan;
%     allCumStrainMap(allCumStrainMap==0) = nan;
%     
%     twinPixelMap(twinPixelMap==0) = nan;
%     twinCumStrainMap(twinCumStrainMap==0) = nan;
%     
%     nonTwinPixelMap(nonTwinPixelMap==0) = nan;    
%     nonTwinCumStrainMap(nonTwinCumStrainMap==0) = nan;

    
    
    try
        close(hWaitbar);
    end
    warning('on','all');
    
    save([summaryPath,'summary_tnont_iE=',num2str(iE),'.mat'],'s_x','s_y','s_euler','s_tp','s_ts','s_ap','s_as','s_fp','s_fs');
end


%% whatever to plot
close all; 
for iE = iE_start:iE_stop
    load([summaryPath,'summary_tnont_iE=',num2str(iE),'.mat']);
    figure('Position',[680   600   750   330]);
    hold on;
    axis equal;axis([0 1.2 0 0.6]);
    axis off;
    theta_d = [0:1:30];  % just the paremeter 'theta'
    border_x = [0 cosd(theta_d) 0];
    border_y = [0 sind(theta_d) 0];
    plot(border_x, border_y, 'k-');
    
%     aData = [s_x, s_y, s_ap, s_fp]; % marker size = num of pixels in grian; color = fraction of pixel by twin 
    aData = [s_x, s_y, s_ap, s_fs]; % marker size = num of pixels in grian; color = fraction of strain by twin 
    
    aData = sortrows(aData,3,'descend');
    ind = (aData(:,3)>0);
    aData = aData(ind,:);
    
    xData = aData(:,1);
    yData = aData(:,2);
    sData = aData(:,3);
    sData = sData/3000;
    cData = aData(:,4);
    
    
    scatter(xData, yData, sData, cData, 'filled');
    
    text(-0.12,0.12,'0 0 0 1','fontsize',18)
    text(1.02,0.02,'2 -1 -1 0','fontsize',18)
    text(0.87,0.52,'1 0 -1 0','fontsize',18)
    colorbar; caxis([0 0.8]);
    set(gca,'fontsize',18);
    hold off;
    output_name = ['cd_ipf nPixels vs twinFraction iE=',num2str(iE)];
    print([summaryPath,output_name,'.tiff'],'-dtiff');
end

%% This can be explained by the SF of extension twins
eu = randi(360,1000,3);
calculate_pos_on_IPF(eu,'Mg','plotTF',1,'ss_of_interest',[19:24],'stressTensor',[-1 0 0; 0 0 0; 0 0 0]);
axis off;
%%
% figure;
% hold on;
% axis square;axis off;
% axis([0 1.2 0 1.2]);
% theta_d = [0:1:30];  % just the paremeter 'theta'
% border_x = [0 cosd(theta_d) 0];
% border_y = [0 sind(theta_d) 0];
% plot(border_x, border_y, 'k-');
% 
% sData = twinPixelMap(:)./allPixelMap(:);
% sData = 300*sData/nanmax(sData(:));
% 
% cData = abs(twinCumStrainMap(:)./ allCumStrainMap(:));
% scatter(xMap(:),yMap(:),sData,cData,'filled')
% 
% text(-0.12,0.12,'0 0 0 1','fontsize',18)
% text(1.02,0.02,'2 -1 -1 0','fontsize',18)
% text(0.87,0.52,'1 0 -1 0','fontsize',18)
% colorbar;
% hold off;
%%
figure;
histogram(eG*10+twinTF);

[~,~,ir] = unique(eG*10+twinTF);
exxP(isnan(exxP)) = 0;
exx21 = accumarray(ir,exxP);












