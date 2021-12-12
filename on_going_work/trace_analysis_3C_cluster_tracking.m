
% chenzhe, 2018-09-17
%
% Similar to previous code. But this is run immediately after performing clustering.  

clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

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
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
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

%% (0) Load cluster number data at all stops

cluster_number_maps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0
cluster_number_maps_cleaned = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load(fullfile(saveDataPath,fName_c2t_result),'clusterNumMap','stru','clusterNumMapCleaned');
    cluster_number_maps{iE} = clusterNumMap;
    cluster_number_maps_cleaned{iE} = clusterNumMapCleaned;
    
    % initialize/zero fields
    for iS =1:length(stru)
        stru(iS).preCluster = zeros(size(stru(iS).cLabel));
        stru(iS).postCluster = zeros(size(stru(iS).cLabel));
        
        stru(iS).gVol = 0;
        stru(iS).volEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).volEvoCleaned = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbMax = zeros(size(stru(iS).cLabel));
    end
    struCell{iE} = stru;
end

%% (1.1) match clusters and track volume evolution, fields related: cVol, cVolCleaned, preCluster, postCluster
for iE = iE_start:iE_stop-1
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(struA)
        % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        
        ID_current = struA(iS).gID;
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        cMapA = cluster_number_maps{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapP = cluster_number_maps{iE+1}(indR_min:indR_max, indC_min:indC_max);
        cMapA_cleaned = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapP_cleaned = cluster_number_maps_cleaned{iE+1}(indR_min:indR_max, indC_min:indC_max);
        
        cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
        cMapP(ID_local~=ID_current) = 0;
        cMapA_cleaned(ID_local~=ID_current) = 0;
        cMapP_cleaned(ID_local~=ID_current) = 0;
        
        struA(iS).gVol = sum(ID_local(:)==ID_current);
        struP(iS).gVol = sum(ID_local(:)==ID_current);
        % match clusters, and find the one with area grown.
        cOverlap = [];
        for ii=1:length(struA(iS).cLabel)
            cNumA = struA(iS).cLabel(ii);
            for jj=1:length(struP(iS).cLabel)
                cNumP = struP(iS).cLabel(jj);
                cOverlap(ii,jj) =sum(sum((cMapA_cleaned==cNumA)&(cMapP_cleaned==cNumP)));   % use cleaned map to calculate how good they overlap !!! But NOT the actual cluster size ------------------
            end
        end
        
        % this is more precise:
        volA_noClean = struA(iS).cVol;
        volA_cleaned = struA(iS).cVolCleaned;
        
        % vrFwd = volP./volA;
        overlapPctA = cOverlap./volA_cleaned;
        % 2018-03-09, I think the reason for the divider to use 'cleaned' is that, the pre-cluster might be real matching, just noisy...
        % If use ./volA_noClean, maybe the 15% threshold can block too many real matched pre-clusters
        
        % do an additional clean up of cOverlap, ignore those which only have <15% overlap with post-cluster
        cOverlap(overlapPctA<0.15) = 0;
        
        [cFrom,cTo] = hungarian_assign(max(cOverlap(:))-cOverlap);
        link = false(size(cOverlap));
        for ii = 1:length(cFrom)
            if (cFrom(ii)>0)&&(cTo(ii)>0)
                link(cFrom(ii),cTo(ii)) = true;
            end
        end
        % twinLikely = link & (overlapPctA>0.5);
        for ii = 1:size(link,1)
            for jj=1:size(link,2)
                if(link(ii,jj))
                    struA(iS).postCluster(ii) = jj;
                    struP(iS).preCluster(jj) = ii;
                end
            end
        end
        
        waitbar(iS/length(stru), hWaitbar);
        
    end
    
    try
        close(hWaitbar);
    end
    warning('on','all');
    
    struCell{iE} = struA;   % pre
    struCell{iE+1} = struP; % post
    
end

% (after 1.1) update to Save. Can disable to temporarily to prevent data loss.
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save(fullfile(saveDataPath,fName_c2t_result),'stru','-append');
end


%% (1.2) after tracking all strain level pairs, get volume history in the whole strain levels, and analyze. Related fields: volEvo, volEvoCleaned, cvInc, (added later: cvIncAfter, cVolGrowthRatio)
%  First, re-load stru in all strain levels
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load(fullfile(saveDataPath,fName_c2t_result),'stru');
    
    % initialize/zero related fields
    for iS =1:length(stru)
        stru(iS).volEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).volEvoCleaned = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbMax = zeros(size(stru(iS).cLabel));
        stru(iS).cVolGrowthRatio = zeros(length(stru(iS).cLabel),length(STOP)-1);
    end
    struCell{iE} = stru;
end


for iE = iE_start:iE_stop
    % stru = struCell{iE};
    for iS =1:length(struCell{iE})
        
        for iCluster = 1:length(struCell{iE}(iS).cLabel)
            
            % [1st loop] find the volume of this area in different strain levels
            
            [iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE, iCluster);
            vol = zeros(1,length(iE_list));
            vol_cleaned = zeros(1,length(iE_list));
            %             tProbEvo = zeros(1,length(iE_list));
            for ii = 1:length(iE_list)
                vol(ii) = struCell{iE_list(ii)}(iS).cVol(iC_list(ii));    % record size of the current cluster in the current iE
                vol_cleaned(ii) = struCell{iE_list(ii)}(iS).cVolCleaned(iC_list(ii));    % record size of the current cluster overlaid with the post cluster, and cleaned
                %                 tProbEvo(ii) = struCell{iE_list(ii)}(iS).tProb(iC_list(ii));
            end
            
            % Method-(3), cVolGrowthRatio = (end_size - smallest_size_up_to_now)/grain_size
            vol_valid = vol_cleaned;
            vol_valid(vol_valid<0) = 0;     % clusters too small, vol/vol_cleaned was made as negative. (Artificial treatment)
            
            cVolGrowthRatio = zeros(size(vol_valid));
            for ii=1:length(vol_valid)
                cVolGrowthRatio(ii) = (vol_valid(end) - min(vol_valid(1:ii)))/struCell{iE}(iS).gVol;
            end
            
            % [2nd loop] copy, consider all iE_to_assign values (should compare with iE, and assign accordingly and differently)
            % (a) if it is related to history, assign to all iE_to_assign values
            % (b) for values depend on current iE and later iEs, use if(iE_to_assign>=iE)
            if length(iE_list) >= 1
                for ii = 1:length(iE_list)
                    iE_to_assign = iE_list(ii);
                    iC_to_assign = iC_list(ii);
                    
                    struCell{iE_to_assign}(iS).volEvo(iC_to_assign,iE_list) = vol;    % record size history of the current cluster in the current iE
                    struCell{iE_to_assign}(iS).volEvoCleaned(iC_to_assign,iE_list) = vol_cleaned;
                    %                     struCell{iE_to_assign}(iS).tProbEvo(iC_to_assign,iE_list) = tProbEvo;
                    
                    struCell{iE_to_assign}(iS).cVolGrowthRatio(iC_to_assign,iE_list) = cVolGrowthRatio;
                    
                    % The search is form iE_start to iE_to_assign !!!
                    %                     struCell{iE_to_assign}(iS).tProbMax(iC_to_assign) = max(tProbEvo(1:find(iE_list==iE_to_assign)));
                    
                end
                
            end
        end
    end
end

% (after 1.2) update to save. Can disable to temporarily to prevent data loss.
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save(fullfile(saveDataPath,fName_c2t_result),'stru','-append');
end


