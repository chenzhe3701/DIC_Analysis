
% chenzhe, 2018-02-27
% based on already labeled data, for twin vs no-twin, look at the criterion
% interested.
%
% (1) iE = xx
% Map = ()
% for each grain, each cluster
%  if twin/notwin, do intereted calculation, add to plot of interested
% Final plot.



% chenzhe, 2018-02-26 based on modified/cleaned cluster result, put cluster
% images into twin and notwin folder.
% Can use as ground truth, but in fact, it still contain mis-identified
% clusters, because visualization can miss a lot of things ...

clear;
addChenFunction;
%%
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

allImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','Choose/make a parent path for output images/files. If have a new task, then make a new one.'),'\'];
mkdir(allImgPath,'twin');
mkdir(allImgPath,'notwin');

img_size = 227; % 227 for alexnet, 224 for vgg, googlenet

useStrain = 0;
makeNewMap = 0;
makeFigure = 0;
summarizeTruth = 1;

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

%% select iE to analyze

% for iE = iE_start:iE_stop
    
    % load data for this iE
    warning('off','all');
    iE = 2;
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
        qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
        qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
        qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
        ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
        exx(ind_outlier) = nan;
        exy(ind_outlier) = nan;
        eyy(ind_outlier) = nan;
    end
    
    
    % use this as the ground truth. Note that the field 'tProb' and the map 'cnnTwinMap' might need to be processed by cnn.
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result]);

    if makeNewMap
        moi = zeros(size(exx)); % new map of interest
    end
    if summarizeTruth
        varTwin = [];
        varNotwin = [];
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
        
        % load previously-made clusterNumMap for this grain
        if makeNewMap
            clusterNumMapLocal = clusterNumMap(indR_min:indR_max, indC_min:indC_max);
            clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
            moiLocal = zeros(size(clusterNumMapLocal));
        end
        
        
        % [code here] to make a figure of something
        if makeFigure
            nCluster = length(stru(iS).cLabel);
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                outputName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
                
                f = figure;
                set(f,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')','ResizeFcn','set(gcf,''visible'',''on'')');
                % plot on figure
                
                
                % write img to parent folder, and then to each labeld folder
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>=0)
                    saveas(f, fullfile(allImgPath,'twin',outputName));
                else
                    saveas(f, fullfile(allImgPath,'notwin',outputName));
                end
            end
        end
        
        
        % [code here] to summarize something for each cluster
        if summarizeTruth
            nCluster = length(stru(iS).cLabel);
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                
                % Apply certain criterion to find the best matching twin system
                pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
                % pdistCS(stru(iS).tSF < sfCF) = nan;             % [criterion-2] SF must > 0.35
                [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
                m_dist = pdistCS(ind_t);
                
                if useStrain
                    % calculate the within-cluster-distance: a vector d
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    moiLocal(indClusterLocal) = 1;
                    
                    data = [exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:))];
                    d = pdist2(data, stru(iS).cCen(iCluster,:));
                    
                    vd = nanmean(d);     % (d) avg within cluster dist
                    ve = nanstd(d);
                end
                
                
                % summarize the varialbe_of_interest, Catogory based on Ground-Truth result, but variable can be calculated in different ways
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>=0)
                    va = m_dist;                 % (a) RSS, disSimilarity, distance
                    va = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain(ind_t,:), 'minkowski', 0.5);
                    
                    vb = stru(iS).tSF(ind_t);    % (b) Schmid factor
                    vc = stru(iS).tProb(iCluster);   % (c) cnnTwinProb
                    
                    vd = stru(iS).vrFwd(iCluster);
                    ve = stru(iS).vrBwd(iCluster);
                    varTwin = [varTwin; ID_current, iCluster, va, vb, vc, vd, ve];
                else
                    va = m_dist;
                    va = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain(ind_t,:), 'minkowski', 0.5);
                    
                    vb = stru(iS).tSF(ind_t);
                    vc = stru(iS).tProb(iCluster);
                    
                    vd = stru(iS).vrFwd(iCluster);
                    ve = stru(iS).vrBwd(iCluster);
                    varNotwin = [varNotwin; ID_current, iCluster, va, vb, vc, vd, ve];
                end
            end
        end
        
        
        % [code here] to make a new map of something
        if makeNewMap
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                
                indClusterLocal = (clusterNumMapLocal==cNum);
                moiLocal(indClusterLocal) = 1;
            end
            
            % copy identified twin system number to twinMap
            moi(indR_min:indR_max, indC_min:indC_max) = moi(indR_min:indR_max, indC_min:indC_max) + moiLocal;
        end
        
        
        
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    try
        close(hWaitbar);
    end
    warning('on','all');
% end


%% whatever to plot
figure;
histogram(varTwin(:,5));
figure;
histogram(varNotwin(:,5));

%%
figure; hold on;
plot3(varTwin(:,3),varTwin(:,4),varTwin(:,6),'ro');
plot3(varNotwin(:,3),varNotwin(:,4),varNotwin(:,6),'bx');
xlabel('RSS');ylabel('SF');zlabel('vInc');

%%
figure; hold on;
plot(varTwin(:,3),varTwin(:,4),'ro');
plot(varNotwin(:,3),varNotwin(:,4),'bx');
xlabel('RSS');ylabel('SF');
%%
fplot(@(x) 3*x+0.35,[0 0.05],'k');

