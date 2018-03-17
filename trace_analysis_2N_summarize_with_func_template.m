
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

% clear;
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

summaryImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','Choose/make a parent path for output images/files. If have a new task, then make a new one.'),'\'];
mkdir(summaryImgPath,'Twin');
mkdir(summaryImgPath,'TwinEnabled');
mkdir(summaryImgPath,'TwinDisabled');
mkdir(summaryImgPath,'Notwin');

twinImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose the parent path for all the source twin images.'),'\'];
mkdir(summaryImgPath,'Classified_twin');
mkdir(summaryImgPath,'Classified_notwin');

img_size = 227; % 227 for alexnet, 224 for vgg, googlenet

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


%% select iE to analyze

useStrain = 0;
makeNewMap = 0;
makeFigure = 0;
summarizeTruth = 1;

for iE = 3%iE_start:iE_stop
    
    % load data for this iE
    warning('off','all');
    % iE = 4;
    
    % use this as the ground truth. Note that the field 'tProb' and the map 'cnnTwinMap' might need to be processed by cnn.
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result], 'stru','clusterNumMap');
    
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
    
    
    if makeNewMap
        map = zeros(size(exx)); % new map of interest
        map2 = zeros(size(exx)); % new map of interest
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
                
%                 % move classified twin images
%                 if stru(iS).tProb(iCluster)>0.5
%                     copyfile(fullfile(twinImgPath,outputName),fullfile(summaryImgPath,'Classified_twin',outputName));
%                 else
%                     copyfile(fullfile(twinImgPath,outputName),fullfile(summaryImgPath,'Classified_notwin',outputName));
%                 end
                
                
                iE_list = iE;
                iC_list = iCluster;
                while 0 ~= struCell{iE_list(1)}(iS).preCluster(iC_list(1))
                    iC_list = [struCell{iE_list(1)}(iS).preCluster(iC_list(1)), iC_list];
                    iE_list = [iE_list(1)-1, iE_list];
                end
                while 0 ~= struCell{iE_list(end)}(iS).postCluster(iC_list(end))
                    iC_list = [iC_list,struCell{iE_list(end)}(iS).postCluster(iC_list(end))];
                    iE_list = [iE_list, iE_list(end)+1];
                end
                
                % write img to parent folder, and then to each labeld folder
                f = figure;
                set(f,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')','ResizeFcn','set(gcf,''visible'',''on'')');
                % plot on figure
                hold on;
                try
                    plot(iE_list, stru(iS).volEvo(iCluster,iE_list), '-ro');
                    plot(iE_list, stru(iS).volEvoCleaned(iCluster,iE_list), '-bx');
                    plot(iE_list, max(stru(iS).volEvo(iCluster,iE_list)) * stru(iS).tProbEvo(iCluster,iE_list), '-gd');
                end
                
                if (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)>0)
                    % save and close
                    saveas(f, fullfile(summaryImgPath,'TwinEnabled',outputName));
                    
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)==0)
                    % save and close
                    saveas(f, fullfile(summaryImgPath,'Twin',outputName));
                    
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)<0)
                    % save and close
                    saveas(f, fullfile(summaryImgPath,'TwinDisabled',outputName));
                else
                    % save and close
                    saveas(f, fullfile(summaryImgPath,'Notwin',outputName));
                end
                close(f);
            end
        end
        
        
        % [code here] to summarize something for each cluster
        if summarizeTruth
            nCluster = length(stru(iS).cLabel);
            for iCluster = 1:nCluster
                cNum = stru(iS).cLabel(iCluster);
                
                % Apply certain criterion to find the best matching twin system
                pdistCS = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain, 'minkowski', 2);       % pair distance between cluster centroids and twinSystem strain components. Non-candidate twinSys lead to nan.
                pdistCS(stru(iS).tSF < 0) = nan;             % [criterion-2] SF must > 0.35
                [m_dist, ind_t] = nanmin(pdistCS,[],2);         % [criterion-3] choose the smallest distanced twinSystem -- [minVal, ind], ind is the corresponding twin system number
                m_dist = pdistCS(ind_t);
                
                
                if useStrain
                    % calculate the within-cluster-distance: a vector d
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    % mapLocal(indClusterLocal);
                    
                    data = [exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:)),exx_local(indClusterLocal(:))];
                    d = pdist2(data, stru(iS).cCen(iCluster,:));
                    
                    vd = nanmean(d);     % (d) avg within cluster dist
                    ve = nanstd(d);
                end
                
                
                % summarize something [varTwin, varNotwin], Catogory based on Ground-Truth result, but variable can be calculated in different ways
                v = [];
                v(1) = ID_current;
                v(2) = iCluster;
                v(3) = m_dist;                 % (a) RSS, disSimilarity, distance
                v(4) = stru(iS).tSF(ind_t);    % (b) Schmid factor
                v(5) = stru(iS).cvInc(iCluster);
                v(6) = stru(iS).tProbMax(iCluster);   % (d) cnnTwinProb
                v(7) = stru(iS).tProb(iCluster);
                
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
                    
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)==0)  % natural twin
                    cNum = stru(iS).cLabel(iCluster);
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    mapLocal(indClusterLocal) = 2;
                elseif (stru(iS).c2t(iCluster)>nss) && (stru(iS).cEnable(iCluster)<0)  % considered as twin, but disabled
                    cNum = stru(iS).cLabel(iCluster);
                    indClusterLocal = (clusterNumMapLocal==cNum);
                    mapLocal(indClusterLocal) = 3;
                end
            end
            
            % copy identified twin system number to twinMap
            map(indR_min:indR_max, indC_min:indC_max) = map(indR_min:indR_max, indC_min:indC_max) + mapLocal;
            %             map2(indR_min:indR_max, indC_min:indC_max) = map2(indR_min:indR_max, indC_min:indC_max) + mapLocal2;
        end
        
        
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    if summarizeTruth
        vAll = [vAll;varTwin;varTwinEnabled;varTwinDisabled;varNotwin];
        vAll(sum(isnan(vAll),2)>0,:)=[];
    end
    
    
    try
        close(hWaitbar);
    end
    warning('on','all');
end


%% whatever to plot
figure; hold on;

n=1;
plot3(varNotwin(:,3),       0.5-varNotwin(:,4)-varNotwin(:,5).*varNotwin(:,6)*n,         varNotwin(:,5).*varNotwin(:,6)*n,           'b.');  % bx

plot3(varTwin(:,3),         0.5-varTwin(:,4)-varTwin(:,5).*varTwin(:,6)*n,           varTwin(:,5).*varTwin(:,6)*n,               'r.','markersize',12);  % ro
plot3(varTwinEnabled(:,3),  0.5-varTwinEnabled(:,4)-varTwinEnabled(:,5).*varTwinEnabled(:,6)*n,    varTwinEnabled(:,5).*varTwinEnabled(:,6)*n, 'm.','markersize',16);  % mo

try
    plot3(varTwinDisabled(:,3),0.5-varTwinDisabled(:,4)-varTwinDisabled(:,5).*varTwinDisabled(:,6)*n,varTwinDisabled(:,5).*varTwinDisabled(:,6)*n,'k.','markersize',12); % co
end

xlabel('RSS');ylabel('SF');zlabel('vol and tProb');


%% a useful way to plot, compare the effect of cvInc*tProbMax, vs individual of them alone  
figure; hold on;
n=1;
opt = 1;
if 1==opt
plot3(varNotwin(:,3),       0.5-varNotwin(:,4),         varNotwin(:,5).*varNotwin(:,6)*n,           'b.');  % bx
plot3(varTwin(:,3),         0.5-varTwin(:,4),           varTwin(:,5).*varTwin(:,6)*n,               'r.','markersize',12);  % ro
plot3(varTwinEnabled(:,3),  0.5-varTwinEnabled(:,4),    varTwinEnabled(:,5).*varTwinEnabled(:,6)*n, 'm.','markersize',16);  % mo
try
    plot3(varTwinDisabled(:,3),0.5-varTwinDisabled(:,4),varTwinDisabled(:,5).*varTwinDisabled(:,6)*n,'k.','markersize',12); % co
end
elseif 1==opt
plot3(varNotwin(:,3),       0.5-varNotwin(:,4),         varNotwin(:,5)*n,           'b.');  % bx
plot3(varTwin(:,3),         0.5-varTwin(:,4),           varTwin(:,5)*n,               'r.','markersize',12);  % ro
plot3(varTwinEnabled(:,3),  0.5-varTwinEnabled(:,4),    varTwinEnabled(:,5)*n, 'm.','markersize',16);  % mo
try
    plot3(varTwinDisabled(:,3),0.5-varTwinDisabled(:,4),varTwinDisabled(:,5)*n,'k.','markersize',12); % co
end
elseif 1==opt
plot3(varNotwin(:,3),       0.5-varNotwin(:,4),         varNotwin(:,6)*n,           'b.');  % bx
plot3(varTwin(:,3),         0.5-varTwin(:,4),           varTwin(:,6)*n,               'r.','markersize',12);  % ro
plot3(varTwinEnabled(:,3),  0.5-varTwinEnabled(:,4),    varTwinEnabled(:,6)*n, 'm.','markersize',16);  % mo
try
    plot3(varTwinDisabled(:,3),0.5-varTwinDisabled(:,4),varTwinDisabled(:,6)*n,'k.','markersize',12); % co
end
end

xlabel('RSS');ylabel('SF');zlabel('vol and or tProb');


% set(gca,'xlim',[0 0.04]);
% 
% set(gca,'ylim',[0 0.5]);
% 
% set(gca,'zlim',[0.01 1.7]);
%%
n = 1;
figure;
plot3(vAll(:,3), 0.5-vAll(:,4)-vAll(:,5).*vAll(:,6)*n, vAll(:,5).*vAll(:,6), 'b.');  % ro
%%
figure; hold on;

plot(varNotwin(:,3),varNotwin(:,4),'bx');

plot(varTwin(:,3),varTwin(:,4),'ro');
plot(varTwinEnabled(:,3),varTwinEnabled(:,4),'mo');

try
    plot(varTwinDisabled(:,3),varTwinDisabled(:,4),'ko');
end
xlabel('RSS');ylabel('SF');
%%
fplot(@(x) 3*x+0.35,[0 0.05],'k');

%%
figure;
histogram(varTwin(:,5));
figure;
histogram(varNotwin(:,5));

%%
myplot(X,Y,map,boundaryTFB); caxis([0,2])
myplot(X,Y,map2,boundaryTFB); caxis([0,2])








