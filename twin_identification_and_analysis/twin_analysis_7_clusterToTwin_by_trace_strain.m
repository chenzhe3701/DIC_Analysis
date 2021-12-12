
% chenzhe, 2018-09-17
%
% This code combines traceDirection analysis and traceStrain analysis to
% identify active twin/slip system.
%
% Run this code after cluster tracking, and finding gb clusters, etc.
%
% chenzhe, 2018-10-20.
% A temporarily ok version.

clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
% dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
% dicFiles = dir([dicPath,'\*.mat']);
% dicFiles = struct2cell(dicFiles);
% dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_2021_09','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
end
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f2 = '_';

debugTF = 0;

%% (0) load data
% store all the clusterNumMap s, omit stop-0
% cluster_number_maps = cell(1,length(STOP)-1);    
clusterNumberMapCell = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load(fullfile(saveDataPath,fName_c2t_result),'stru','clusterNumMap','clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
    twinMapCell{iE} = zeros(size(clusterNumMapCleaned));
    sfMapCell{iE} = zeros(size(clusterNumMapCleaned));
    cToGbDistMapCell{iE} = zeros(size(clusterNumMapCleaned));
    % initialize/zero related fields
    for iS =1:length(stru)
        stru(iS).cActiveSS = zeros(length(stru(iS).cLabel), length(stru(iS).tLabel));
    end
    % try to remove some fields, if needed
    try
        stru = rmfield(stru,{'tR2'});
    end
    struCell{iE} = stru;
end

%% a script to summarize typical number of peaks. added to revise paper part-1, 2020-07-07.
if 0
    trace_analysis_3D_aux_summarize_npeaks();
    save('peakInfo.mat','peakInfo')
    figure;
    histogram(peakInfo(:,5),0:50:1500);
    xlabel('Number of Peaks');
    ylabel('Counts');
    set(gca,'fontsize',18)
end

%% (1) analyze
warning('off','all');

% tunning parameters
p.diffStrain_cr = 0.05;   % diff in strain need to < 0.07;
p.rEffStrain = [0.2, 1.8];    % ratio of cluster effective strain and twin effective strain should be within this range   
p.strainRank_cr = 0.3;  % rank need to >=1 (in 0-based) 
p.SF_th = 0.2;        % twin SF need to be larger than this  (or 0.1)
p.pctTotalPeaks_bwd_cr = 0.2; % when going backward, need at least these pct of total number of peaks to be considered as active trace  (or 0.05, 0.15)  
% note that this is 'th_1' input outside of this function  
p.pctTotalPeaks_fwd_cr = 0.3;
p.clusterToGbPct_th = 0.2; % 95% quantile of cluster point's distance to gb. if < pct * gDia, considered as too close, and not identify as twin.

pctTotalPeaks_fwd_cr = p.pctTotalPeaks_fwd_cr; 
                    
stru = struCell{iE_start};
for iS = 1:length(stru)
    %     iS = find(arrayfun(@(x) x.gID == 179,stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace]
    %     iS = find(gIDwithTrace == 296); % for debugging.
    close all;
    ID_current = stru(iS).gID;
    
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
    
    % (2) Here, we want to do analysis on clustered maps, to decide the active slip/twin system.
    % For examples in twinning analysis, we previously performed cluster analysis and 'identified'/'confirmed' twin clusters.
    % More generally, we might need to first do a rough clustering based on strain map, then perform trace analysis, to decide which are the exist slip/twin systems.
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]); 
    [nR,nC] = size(ID);
    indC_min = find(sum(ind_local, 1), 1, 'first'); 
    indC_max = find(sum(ind_local, 1), 1, 'last');  
    indR_min = find(sum(ind_local, 2), 1, 'first'); 
    indR_max = find(sum(ind_local, 2), 1, 'last');  
    % increase area size by 1 pixel. 2021-09-29
    indC_min = max(indC_min-1, 1);
    indC_max = min(indC_max+1, nC);
    indR_min = max(indR_min-1, 1);
    indR_max = min(indR_max+1, nR);
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
    
    % initialize for each grain (iS)
    for iE = iE_start:iE_stop
        twinMapLocal{iE} = zeros(size(ID_local));
        sfMapLocal{iE} = zeros(size(ID_local));
        cToGbDistMapLocal{iE} = zeros(size(ID_local));
    end
    twinMapCell_cluster = []; % for cluster
    sfMapCell_cluster = [];
    
    % for each iE_entry (the entry point for analysis)
    for iE_entry = iE_start:iE_stop
        % for each iC_outer
        for iC_entry = 1:length(struCell{iE_entry}(iS).cLabel)
            
            % We need to analyze this cluster [iC_outer] at the strain level [iE_outer], but this will need the information from the tracked [iE_list] and [iC_list].
            % So, first find the [iE_list, iC_list]
            [iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_entry, iC_entry);
            
            % Analyze all the linked iEs.  So, if iE_list(1)==iE_outer, it means it has not been analyzed before, then do [iE_list(ii),iC_list(ii)] pairs
            if iE_list(1) == iE_entry  
                for iEC = 1:length(iE_list)
                    close all;
                    iE = iE_list(iEC);
                    iC = iC_list(iEC);
                    
                    clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
                    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
                    if debugTF
                        myplot(clusterNumMapL);
                    end
                    
                    ssAllowed = ones(ntwin,1);
                    
%                     [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = ...
%                         label_twin_trace(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
%                         struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,'twin',debugTF, pctTotalPeaks_fwd_cr, 0.3, ssAllowed);
                    [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = ...
                        label_twin_trace_with_stats(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local,uniqueBoundary_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
                        struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceDir,traceND,traceSF,sampleMaterial,'twin',debugTF, pctTotalPeaks_fwd_cr, p, ssAllowed);
                    % each cell contains cells of tMap at an iEs
                    
                end % end of iEC
                
            end
            
        end % end of iC_outer
        
    end % end of iE_outer
    
    % for each strain level, update twinMapLocal{iE} with tMapCell
    for iE = iE_start:iE_stop
        for jj = 1:size(twinMapCell_cluster,2)
            if ~isempty(twinMapCell_cluster{iE,jj})
                twinMapLocal{iE} = twinMapLocal{iE} + twinMapCell_cluster{iE,jj};
                sfMapLocal{iE} = sfMapLocal{iE} + sfMapCell_cluster{iE,jj};
            end
        end
        % double check ================================================================================================================> 2021-09-28
        if any(~ismember(unique(twinMapLocal{iE}(:)),0:6))
            error(' ')
        end
        
        clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
        cToGbDistMapLocal{iE} = zeros(size(ID_local));
        for iC = 1:length(struCell{iE}(iS).cLabel)
            cNum = struCell{iE}(iS).cLabel(iC);
            cToGbDistMapLocal{iE}(clusterNumMapL==cNum) = struCell{iE}(iS).cToGbDist(iC,end);
        end
        
        % update. First clean old map, then add new map.
        map_local = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);  % (1) Cut a squared map from big map
        twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;          % (2) Eliminate this squared region from the big map
        map_local(ID_local == ID_current) = 0;                              % (3) clean the grain area on the cut map
        % ID_local(128,37)
        map_local = map_local + twinMapLocal{iE};                           % (4) update the grain area on the cut map
        
        % double check ===========================================================================================================================> 2021-09-28
        if any(~ismember(unique(map_local(:)),0:6))
            
            error(' ')
        end
        twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;  % (5) Add the modified cut map to big map
        
        map_local = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max); 
        sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;
        map_local(ID_local == ID_current) = 0; 
        map_local = map_local + sfMapLocal{iE}; 
        sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;
        
        map_local = cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;
        map_local(ID_local == ID_current) = 0; 
        map_local = map_local + cToGbDistMapLocal{iE};
        cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;
    end
    disp(iS);
end % end of iS
warning('on','all');

save(fullfile(saveDataPath,[sampleName,'_1_twin_result_by_trace_strain.mat']),'twinMapCell','sfMapCell','cToGbDistMapCell','struCell','p','-v7.3');
%% Need more code to analyze accuracy.  
% Note 2021-09-26 when recovering from disk crash 2021-08: 
% this seems like repeating the analysis of assessing identification parameters for MaterChar tiwn identification paper 2020.
% Some TP/FP/FN/TN analysis.

if 0
    
[truthFile, truthPath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map.mat','select the truth results for twin-grain boundary intersection');
[checkFile, checkPath] = uigetfile('D:\p\m\DIC_Analysis\20191230_1756_twinMaps_K.mat','select the results for twin-grain boundary intersection to check');
d = load(fullfile(truthPath,truthFile),'struCell');
trueStruCell = d.struCell;
load(fullfile(truthPath,truthFile),'trueTwinMapCell');

% In some codes twin systems use number 19:24. 
% but variant map -> trueTwinMap uses 1:6.
% So, check first.
if max(trueTwinMapCell{3}(:))>6
    for iE=2:5
        temp = trueTwinMapCell{iE};
        temp(temp>18) = temp(temp>18)-18;
        trueTwinMapCell{iE} = temp;
    end
end

load(fullfile(checkPath,checkFile),'struCell');
load(fullfile(checkPath,checkFile),'twinMapCell');
if max(twinMapCell{3}(:))>6
    for iE=2:5
        temp = twinMapCell{iE};
        temp(temp>18) = temp(temp>18)-18;
        twinMapCell{iE} = temp;
    end
end

%% (1) count summary (grain level)
mt_all = [];
mi_all = [];
conf_mat = [];
for iE = 2:5
    mt = [];
    mi = [];
   for iS = 1:length(struCell{iE})
       numC = length(struCell{iE}(iS).cLabel);
       % number accuracy
       cTrueTwin = trueStruCell{iE}(iS).cTrueTwin;
       mt = [mt;cTrueTwin];
       
       cActiveSS = struCell{iE}(iS).cActiveSS;
       mi = [mi;cActiveSS];
   end
   mt_all = [mt_all; mt];
   mi_all = [mi_all; mi];
   
%    tp = sum((mt(:)==1)&(mi(:)==1));
%    fp = sum((mt(:)==0)&(mi(:)==1));
%    fn = sum((mt(:)==1)&(mi(:)==0));
%    tn = sum((mt(:)==0)&(mi(:)==0));
   
   % method new - 1
   tp = sum(sum((mt>0)&(mi>0)&(mt==mi)));   % should be the same as above  
   fp = sum(sum((mt==0)&(mi>0)));   % should be the same as above
   fn = sum(sum((mt>0)&(mi==0)));   % should be the same as above  
   tn = sum((sum(mt,2)==0)&(sum(mi,2)==0));     % should be different
   
   % method new - 2. 
   % actually twinned = any(mt,2) 
   % identified the same = ~any(mt-mi,2).   
   % identified different = any(mt-mi,2). 
   tp = sum( any(mt,2) & ~any(mt-mi,2) );
   fp = sum( ~any(mt,2) & any(mt-mi,2) );
   fn = sum( any(mt,2) & any(mt-mi,2) );
   tn = sum( ~any(mt,2) & ~any(mt-mi,2) );
   
   conf_mat = [conf_mat;tp,fp,fn,tn];
end
conf_mat_all = sum(conf_mat,1);

tbl_1 = summarize_confusion(conf_mat);

%% (2) area summary (not cleaned, so a rough summary)
conf_matA = [];
for iE=2:5
    ind = (trueTwinMapCell{iE}==twinMapCell{iE})&(trueTwinMapCell{iE}>0);
    tp = sum(ind(:));
    ind = (trueTwinMapCell{iE}==0)&(twinMapCell{iE}>0);
    fp = sum(ind(:));
    ind = (trueTwinMapCell{iE}>0)&(twinMapCell{iE}==0);
    fn = sum(ind(:));
    ind = (trueTwinMapCell{iE}==twinMapCell{iE})&(trueTwinMapCell{iE}==0);
    tn = sum(ind(:));
    conf_matA = [conf_matA;tp,fp,fn,tn];
end
conf_matA_all = sum(conf_matA,1);

tbl_2 = summarize_confusion(conf_matA);

%% (3) area summary cleaned (pixel level). May need to run codes following this part to clean first.

load(fullfile(checkPath,checkFile),'twinMapCleanedCell');
if max(twinMapCleanedCell{3}(:))>6
    for iE=2:5
        temp = twinMapCleanedCell{iE};
        temp(temp>18) = temp(temp>18)-18;
        twinMapCleanedCell{iE} = temp;
    end
end

conf_matB = [];
confMap = [];
confMapCell = [];
for iE=2:5
    confMap = zeros(size(trueTwinMapCell{iE}));
    
    ind = (trueTwinMapCell{iE}==twinMapCleanedCell{iE})&(trueTwinMapCell{iE}>0);
    confMap(ind) = 3;
    tp = sum(ind(:));
    
    ind = (trueTwinMapCell{iE}==0)&(twinMapCleanedCell{iE}>0);
    confMap(ind) = 2;
    fp = sum(ind(:));
    
    ind = (trueTwinMapCell{iE}>0)&(twinMapCleanedCell{iE}~=trueTwinMapCell{iE});
    confMap(ind) = 1;
    fn = sum(ind(:));
    
    ind = (trueTwinMapCell{iE}==twinMapCleanedCell{iE})&(trueTwinMapCell{iE}==0);
    confMap(ind) = 0;
    tn = sum(ind(:));
    
    confMapCell{iE} = confMap;
    conf_matB = [conf_matB;tp,fp,fn,tn];
end
conf_matB_all = sum(conf_matB,1);

tbl_3 = summarize_confusion(conf_matB);

end

if 0
%% plot confusion map (pixel level summary)
iE_start = 2;
iE_select = 4;
colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];
[f,a,c] = myplot(X, Y, confMapCell{iE_select}, boundaryTFB, 3);
colormap(colorMap);
caxis([-1.5 3.5]);set(c,'limits',[-0.5, 3.5]);
c.Ticks = 0:3;
c.TickLabels={['TN'], ['FN'],['FP'],['TP']};
tp = conf_matB(iE_select-iE_start+1,1);
fp = conf_matB(iE_select-iE_start+1,2);
fn = conf_matB(iE_select-iE_start+1,3);
tn = conf_matB(iE_select-iE_start+1,4);
nPix = tp+fp+fn+tn;
c.TickLabels={['TN: ',[num2str(tn/nPix*100,'%0.2f'),'%']], ['FN: ',[num2str(fn/nPix*100,'%0.2f'),'%']],...
    ['FP: ',[num2str(fp/nPix*100,'%0.2f'),'%']],['TP: ',[num2str(tp/nPix*100,'%0.2f'),'%']]};
% c.TickLabels={['TN: ',num2str(TN)], ['FN: ',num2str(FN)],['FP: ',num2str(FP)],['TP: ',num2str(TP)]};
set(a,'fontsize',18);

%% maximize plot and run this:
script_make_double_axis;
%%
print(fullfile(saveDataPath, [ttl,'.tif']),'-dtiff');
end

%%

if 0
%% cleanup, plot trueTwinMaps, indicate TF/FP/FN/TN values
% select data that want to do clean up
[toCleanFile, toCleanPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the results to do clean up');
load(fullfile(toCleanPath,toCleanFile),'struCell','twinMapCell');

useParallel = 1;
if useParallel    
    for iE = iE_start:iE_stop
        variantMap = twinMapCell{iE};
        
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
                    disp(['ID = ',num2str(ID_current)]);
                end
            end
        end
        
        variantMapCleaned = zeros(size(ID)); % new map of interest
        for ii=1:npool
            variantMapCleaned = variantMapCleaned + partMap{ii};
        end
        
        twinMapCleanedCell{iE} = variantMapCleaned;
    end    
end
save(fullfile(toCleanPath,toCleanFile),'twinMapCleanedCell','-append');


%% 
tp = a(:,1);
fp = a(:,2);
fn = a(:,3);
tn = a(:,4);

PPV = tp./(tp+fp);
TPR = tp./(tp+fn);
ACC = (tp+tn)./(tp+fp+fn+tn);

a(:,5:7) = [ACC, PPV, TPR];



end

