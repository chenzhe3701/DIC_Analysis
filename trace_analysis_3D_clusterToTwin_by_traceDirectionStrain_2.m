
% chenzhe, 2018-09-17
%
% This code combines traceDirection analysis and traceStrain analysis to
% identify active twin/slip system.
%
% Run this code after cluster tracking.

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

%% (0) load data
cluster_number_maps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0
cluster_number_maps_cleaned = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    cluster_number_maps{iE} = clusterNumMap;
    cluster_number_maps_cleaned{iE} = clusterNumMapCleaned;
    twinMap{iE} = zeros(size(clusterNumMapCleaned));
    sfMap{iE} = zeros(size(clusterNumMapCleaned));
    cToGbDistMap{iE} = zeros(size(clusterNumMapCleaned));
    % initialize/zero related fields
    for iS =1:length(stru)
        stru(iS).tR2 = zeros(length(stru(iS).cLabel),length(stru(iS).tLabel));
        stru(iS).cActiveSS = zeros(length(stru(iS).cLabel), length(stru(iS).tLabel));
    end
    struCell{iE} = stru;
end



%% (1) analyze
warning('off','all');
stru = struCell{iE_start};
for iS = 1:length(stru)
    %
    %     iS = find(arrayfun(@(x) x.gID == 378,stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace]
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
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    
    % initialize for each grain (iS)
    for iE = iE_start:iE_stop
        twinMapLocal{iE} = zeros(size(ID_local));
        sfMapLocal{iE} = zeros(size(ID_local));
        cToGbDistMapLocal{iE} = zeros(size(ID_local));
    end
    twinMapCell = [];
    sfMapCell = [];
    r2MapCell = []; % but not used currently
    
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
                    
                    clusterNumMapL = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
                    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
                    if debugTF
                        myplot(clusterNumMapL);
                    end
                    
                    ssAllowed = ones(ntwin,1);
                    [twinMapCell, sfMapCell, r2MapCell, struCell, haveActiveSS] = label_twin_trace(twinMapCell, sfMapCell, r2MapCell, cluster_number_maps_cleaned,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
                        struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,'twin',debugTF, 0.3,0.3,ssAllowed);
                    % each cell contains cells of tMap at an iEs

                end % end of iEC
                
            end
            
        end % end of iC_outer
        
    end % end of iE_outer
    
    % for each strain level, update twinMapLocal{iE} with tMapCell
    for iE = iE_start:iE_stop
        for jj = 1:size(twinMapCell,2)
            if ~isempty(twinMapCell{iE,jj})
                twinMapLocal{iE} = twinMapLocal{iE} + twinMapCell{iE,jj};
                sfMapLocal{iE} = sfMapLocal{iE} + sfMapCell{iE,jj};
            end
        end
        
        clusterNumMapL = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
        cToGbDistMapLocal{iE} = zeros(size(ID_local));
        for iC = 1:length(struCell{iE}(iS).cLabel)
            cNum = struCell{iE}(iS).cLabel(iC);
            cToGbDistMapLocal{iE}(clusterNumMapL==cNum) = struCell{iE}(iS).cToGbDist(iC,end);
        end
        
        % update. First clean old map, then add new map.
        toClean = twinMap{iE}(indR_min:indR_max, indC_min:indC_max);
        toClean(ID_local ~= ID_current) = 0;
        twinMap{iE}(indR_min:indR_max, indC_min:indC_max) = twinMap{iE}(indR_min:indR_max, indC_min:indC_max) + twinMapLocal{iE};
        
        toClean = sfMap{iE}(indR_min:indR_max, indC_min:indC_max);
        toClean(ID_local ~= ID_current) = 0;
        sfMap{iE}(indR_min:indR_max, indC_min:indC_max) = sfMap{iE}(indR_min:indR_max, indC_min:indC_max) + sfMapLocal{iE};
        
        toClean = cToGbDistMap{iE}(indR_min:indR_max, indC_min:indC_max);
        toClean(ID_local ~= ID_current) = 0;
        cToGbDistMap{iE}(indR_min:indR_max, indC_min:indC_max) = cToGbDistMap{iE}(indR_min:indR_max, indC_min:indC_max) + cToGbDistMapLocal{iE};
    end
    disp(iS);
end % end of iS
warning('on','all');

timeStr = datestr(now,'yyyymmdd_HHMM');
save([timeStr,'_twinMaps.mat'],'twinMap','sfMap','cToGbDistMap','struCell','-v7.3');
%%





