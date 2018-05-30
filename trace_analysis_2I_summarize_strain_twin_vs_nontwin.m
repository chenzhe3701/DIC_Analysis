
% chenzhe, 2018-5-26
% based on trace_analysis_2N_().
% reduce redundant codes, just to summarize something of interest.
% (1) The strain (e.g., exx) accommodated by twin vs nontwin, at differetn global strain, using box plot. 
% We dicide not to show this plot, so this code turned out to be historical     
% 
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


%% (0) initialize something
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    struCell{iE} = stru;
end

for iE = iE_start:iE_stop
    
    % load data for this iE
    warning('off','all');
    
    % use this as the ground truth. Note that the field 'tProb' and the map 'cnnTwinMap' might need to be processed by cnn.
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result], 'stru','clusterNumMapCleaned','trueTwinMap');
    
    if useStrain
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
end
%% select iE to analyze

useStrain = 1;
makeNewMap = 0;
makeFigure = 0;
summarizeTruth = 0;

exxP = []; % pixel-wise
exyP = []; % pixel-wise
eyyP = []; % pixel-wise
eeffP = []; % pixel-wise

infoC = []; % cluster-wise
infoG = []; % grain-wise
twinTF = [];
eG = [];
for iE = iE_start:iE_stop
    
    % load data for this iE
    warning('off','all');
    
    % use this as the ground truth. Note that the field 'tProb' and the map 'cnnTwinMap' might need to be processed by cnn.
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result], 'stru','clusterNumMapCleaned','trueTwinMap');
    
    if useStrain
%         strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
%         load(strainFile,'exx','exy','eyy','sigma');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------
%         clear('exy_corrected');
%         load(strainFile,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
%         
%         if exist('exy_corrected','var')&&(1==exy_corrected)
%             disp('================= exy already corrected ! ========================');
%             exy_corrected = 1;
%         else
%             disp('================= exy being corrected here ! =======================');
%             exy = -exy;
%             exy_corrected = 1;
%         end
%         % remove bad data points
%         exx(sigma==-1) = nan;
%         exy(sigma==-1) = nan;
%         eyy(sigma==-1) = nan;
%         eeff = calculate_effective_strain(exx,exy,eyy);
%         qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
%         qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
%         qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
%         qt_eeff = quantile(eeff(:),[0.0013,0.9987]); qt_eeff(1)=min(-1,qt_eeff(1)); qt_eeff(2)=max(1,qt_eeff(2));
%         ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
%         exx(ind_outlier) = nan;
%         exy(ind_outlier) = nan;
%         eyy(ind_outlier) = nan;
%         eeff(ind_outlier) = nan;
        
        exx = EXX{iE};
        exy = EXY{iE};
        eyy = EYY{iE};
        eeff = EEFF{iE};
    end
    
    % down-sampling so that memory can handle
    rr = 25;
    exx_r = exx(1:rr:end,1:rr:end);
    exy_r = exy(1:rr:end,1:rr:end);
    eyy_r = eyy(1:rr:end,1:rr:end);
    eeff_r = eeff(1:rr:end,1:rr:end);
    trueTwinMap_r = trueTwinMap(1:rr:end,1:rr:end);
   
    % (1) Look at exx
    % twin
    ind = (trueTwinMap_r>18);
    vec = exx_r(ind);
    exxP = [exxP;vec];
    vec = exy_r(ind);
    exyP = [exyP;vec];
    vec = eyy_r(ind);
    eyyP = [eyyP;vec];
    vec = eeff_r(ind);
    eeffP = [eeffP;vec];
    twinTF = [twinTF; 1*ones(length(vec),1)];
    eG = [eG; iE*ones(length(vec),1)];
    % nontwin    
    ind = (trueTwinMap_r<18);
    vec = exx_r(ind);
    exxP = [exxP;vec];
    vec = exy_r(ind);
    exyP = [exyP;vec];
    vec = eyy_r(ind);
    eyyP = [eyyP;vec];
    vec = eeff_r(ind);
    eeffP = [eeffP;vec];
    twinTF = [twinTF; 2*ones(length(vec),1)];
    eG = [eG; iE*ones(length(vec),1)];
   
    
    
    if makeNewMap
        map = zeros(size(exx)); % new map of interest
    end
    if summarizeTruth
        varTwin = [];
        varTwinEnabled = [];
        varTwinDisabled = [];
        varNotwin = [];
        varAll = [];
    end
    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(stru)
        % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        % iS = find(gIDwithTrace == 296); % for debugging.
        ID_current = gIDwithTrace(iS);
        s = stru(iS);   % info in stru for this grain
        
        if (useStrain)||(makeNewMap)
            ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
            indC_min = find(sum(ind_local, 1), 1, 'first');
            indC_max = find(sum(ind_local, 1), 1, 'last');
            indR_min = find(sum(ind_local, 2), 1, 'first');
            indR_max = find(sum(ind_local, 2), 1, 'last');
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        end
        
        if useStrain
            exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
            exy_local = exy(indR_min:indR_max, indC_min:indC_max);
            eyy_local = eyy(indR_min:indR_max, indC_min:indC_max);
            boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
            x_local = X(indR_min:indR_max, indC_min:indC_max);
            y_local = Y(indR_min:indR_max, indC_min:indC_max);
            
            
            exx_local(ID_local~=ID_current) = nan;  % make eij anum only for the current grain
            exy_local(ID_local~=ID_current) = nan;
            eyy_local(ID_local~=ID_current) = nan;
            
            % % additional down sampling
            % rr  = 5;
            % exx_local = exx_local(1:rr:end,1:rr:end);
            % exy_local = exy_local(1:rr:end,1:rr:end);
            % eyy_local = eyy_local(1:rr:end,1:rr:end);
            % boundaryTF_local = boundaryTF_local(1:rr:end,1:rr:end);
            % x_local = x_local(1:rr:end,1:rr:end);
            % y_local = y_local(1:rr:end,1:rr:end);
            % ID_local = ID_local(1:rr:end,1:rr:end);
            % %
            
            % find vectors for cluster, using ind
            ind = find((ID_local==ID_current)); %&(~isnan(exx_local)));
            exx_t = exx_local(ind);
            exy_t = exy_local(ind);
            eyy_t = eyy_local(ind);
            data_t = [exx_t(:), exy_t(:), eyy_t(:)];
            
            clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
            clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
        end
        
        % [code here] to make a figure of something
        if makeFigure
            nCluster = length(stru(iS).cLabel);
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                outputName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
            end
        end
        
        
        % [code here] to summarize something for each cluster
        if summarizeTruth
            nCluster = length(stru(iS).cLabel);
            for iCluster = 1:nCluster
                idx = find(gID==ID_current);
                cNum = stru(iS).cLabel(iCluster);
                tsNum = stru(iS).ts(iCluster);
                trueTwin = stru(iS).trueTwin(iCluster);
                dis = stru(iS).dis(iCluster);
                sf = stru(iS).sf(iCluster);
                phi_classifier = 1/8*sf + 7/8*min(dis,1) + 15/16;
                Psize = stru(iS).cVolGrowthRatio(iCluster,iE);
                Pcnn = stru(iS).tProbMax(iCluster);
                eta_classifier = Psize * Pcnn;
                cVolCleaned = stru(iS).cVolCleaned(iCluster);
                gVol = stru(iS).gVol;
                euler = [gPhi1(idx),gPhi(idx),gPhi2(idx)];

                
                if useStrain
                    % calculate the within-cluster-distance: a vector d
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    % mapLocal(indClusterLocal);
                    data = [exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:))];
                end
                
                % summarize something [varTwin, varNotwin], Catogory based on Ground-Truth result, but variable can be calculated in different ways
                v = [];
                v(1) = ID_current;
                v(2) = iCluster;

                
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)==0)
                    varTwin = [varTwin; v];
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>0)
                    varTwinEnabled = [varTwinEnabled; v];
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)<0)
                    varTwinDisabled = [varTwinDisabled; v];
                else
                    varNotwin = [varNotwin; v];
                end
            end

        end
        
        
        % [code here] to make a new map of something
        if makeNewMap
            % load previously-made clusterNumMap for this grain
            clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
            clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
            mapLocal = zeros(size(clusterNumMapLocal));
            for iCluster = 1:length(stru(iS).cLabel)
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>0)   % enabled twin
                    cNum = stru(iS).cLabel(iCluster);
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    mapLocal(indClusterLocal) = 1;
                end
            end
            % copy identified twin system number to twinMap
            map(indR_min:indR_max, indC_min:indC_max) = map(indR_min:indR_max, indC_min:indC_max) + mapLocal;
        end
        
        
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    if summarizeTruth
        varAll = [varAll;varTwin;varTwinEnabled;varTwinDisabled;varNotwin];
        varAll(sum(isnan(varAll),2)>0,:)=[];
    end
    
    
    try
        close(hWaitbar);
    end
    warning('on','all');
end


%% whatever to plot
close all;
figure;
myboxplot(exxP, eG*10+twinTF, 'whisker',[0.03,0.97],'colors','rb','Positions',[1,2, 4,5, 7,8, 10,11],'width',0.8);
set(gca,'ylim',[-0.12,0.02],'fontsize',18,'xTick',[1.5, 4.5, 7.5, 10.5],'xTickLabel',{'-0.004','-0.012','-0.023','-0.039'});
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('\fontsize{24}\epsilon_{xx}');

figure;
myboxplot(exyP, eG*10+twinTF, 'whisker',[0.03,0.97],'colors','rb','Positions',[1,2, 4,5, 7,8, 10,11],'width',0.8);
set(gca,'ylim',[-0.05,0.05],'fontsize',18,'xTick',[1.5, 4.5, 7.5, 10.5],'xTickLabel',{'-0.004','-0.012','-0.023','-0.039'});
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('\fontsize{24}\epsilon_{xy}');

figure;
myboxplot(eyyP, eG*10+twinTF, 'whisker',[0.03,0.97],'colors','rb','Positions',[1,2, 4,5, 7,8, 10,11],'width',0.8);
set(gca,'ylim',[-0.02,0.06],'fontsize',18,'xTick',[1.5, 4.5, 7.5, 10.5],'xTickLabel',{'-0.004','-0.012','-0.023','-0.039'});
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('\fontsize{24}\epsilon_{yy}');
figure;

myboxplot(eeffP, eG*10+twinTF, 'whisker',[0.03,0.97],'colors','rb','Positions',[1,2, 4,5, 7,8, 10,11],'width',0.8);
set(gca,'ylim',[0,0.12],'fontsize',18,'xTick',[1.5, 4.5, 7.5, 10.5],'xTickLabel',{'-0.004','-0.012','-0.023','-0.039'});
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('\fontsize{24}\epsilon_{eff}');

%%
figure;
histogram(eG*10+twinTF);

[~,~,ir] = unique(eG*10+twinTF);
exxP(isnan(exxP)) = 0;
exx21 = accumarray(ir,exxP);












