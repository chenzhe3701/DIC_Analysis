% track cluster NumEvolution and VolEvolution
% chenzhe, 2018-03-01


clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));


% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');

%% (0) initialize something
cNumMaps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0
cNumMaps_cleaned = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    cNumMaps{iE} = clusterNumMap;
    cNumMaps_cleaned{iE} = clusterNumMapCleaned;
    % remove some fields
    try
        stru = rmfield(stru,{'vrFwd','vrBwd'});
    end
    % initialize/zero fields
    for iS =1:length(stru)
        stru(iS).cVol = zeros(size(stru(iS).cLabel));
        stru(iS).preCluster = zeros(size(stru(iS).cLabel));
        stru(iS).postCluster = zeros(size(stru(iS).cLabel));
        stru(iS).preVol = zeros(size(stru(iS).cLabel));
        stru(iS).postVol = zeros(size(stru(iS).cLabel));
        stru(iS).preVolCleaned = zeros(size(stru(iS).cLabel));
        stru(iS).postVolCleaned = zeros(size(stru(iS).cLabel));
        stru(iS).cLikeT = zeros(size(stru(iS).cLabel));     % based on pair comparison between a pair of strain levels
        stru(iS).volEvo = zeros(1,length(STOP)-1);
        stru(iS).cvInc = false(size(stru(iS).cLabel));      % based on looking at the volume evolution in all strain levels
    end
    struCell{iE} = stru;
end

%% (1) match clusters and track volume evolution
for iE = iE_start:iE_stop-1
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(struA)
        % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        % iS = find(gIDwithTrace == 296); % for debugging.
        ID_current = gIDwithTrace(iS);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        cMapA = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapP = cNumMaps{iE+1}(indR_min:indR_max, indC_min:indC_max);
        cMapA_cleaned = cNumMaps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapP_cleaned = cNumMaps_cleaned{iE+1}(indR_min:indR_max, indC_min:indC_max);
        
        cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
        cMapP(ID_local~=ID_current) = 0;
        cMapA_cleaned(ID_local~=ID_current) = 0;
        cMapP_cleaned(ID_local~=ID_current) = 0;
        
        % match clusters, and find the one with area grown.
        cOverlap = [];
        for ii=1:length(struA(iS).cLabel)
            cNumA = struA(iS).cLabel(ii);
            for jj=1:length(struP(iS).cLabel)
                cNumP = struP(iS).cLabel(jj);
                struA(iS).cVol(jj) = sum(cMapA(:)==cNumA);  % count the volume for the current cluster
                cOverlap(ii,jj) =sum(sum((cMapA_cleaned==cNumA)&(cMapP_cleaned==cNumP)));   % use cleaned map to calculate !!! ------------------
            end
        end
        volA_cleaned = sum(cOverlap,2);
        volP_cleaned = sum(cOverlap,1);
        % vrFwd = volP./volA;
        overlapPctA = cOverlap./volA_cleaned;
        
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
                    struA(iS).cLikeT(ii) = 1;
                    
                    struA(iS).postCluster(ii) = jj;
                    struP(iS).preCluster(jj) = ii;
                    
                    struA(iS).postVol(ii) = sum(cMapP(:)==jj);
                    struP(iS).preVol(jj) = sum(cMapA(:)==ii);
                    
                    struA(iS).postVolCleaned(ii) = volP_cleaned(jj);
                    struP(iS).preVolCleaned(jj) = volA_cleaned(ii);
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



%% (2) after tracking all strain level pairs, get volume history in the whole strain levels, and analyze
for iE = iE_start:iE_stop
    stru = struCell{iE};
    for iS =1:length(struCell{iE})
        
        for iCluster = 1:length(struCell{iE}(iS).cLabel)
            
            vol = zeros(1,length(STOP)-1);
            
            % [1st loop] find the volume of this area in different strain levels
            iE_current = iE;
            iC_current = iCluster;
            if 0==struCell{iE_current}(iS).preCluster(iC_current)
                vol(iE_current) = struCell{iE_current}(iS).cVol(iC_current);    % record size of the current cluster in the current iE
                % update iC, iE
                iC_current = struCell{iE_current}(iS).postCluster(iC_current);
                iE_current = iE_current + 1;
            else
                iC_current = 0;
                iE_current = iE_current + 1;                
            end
            while iC_current
                vol(iE_current) = struCell{iE_current}(iS).cVol(iC_current);   % update volume to the cluster size in the next iE
                % update iC, iE
                iC_current = struCell{iE_current}(iS).postCluster(iC_current);
                iE_current = iE_current + 1;                
            end
            
            % evaluate if volume seems to be increasing
            vol_valid = vol(iE:iE_current-1);
            if 1==length(vol_valid)
                vInc = false;
            elseif 2==length(vol_valid)
                vInc = diff(vol_valid)>0;
            else
                vol_valid = vol_valid(2:end);
                kb = ([1:length(vol_valid);ones(1,length(vol_valid))]')\(vol_valid');
                vInc = kb(1)>0;
            end
            
            % [2nd loop] copy the volume to all corresponding clusters at different strain levels
            iE_current = iE;
            iC_current = iCluster;
            if 0==struCell{iE_current}(iS).preCluster(iC_current)
                struCell{iE_current}(iS).volEvo(iC_current,:) = vol;    % record size history of the current cluster in the current iE
                struCell{iE_current}(iS).cvInc(iC_current) = vInc;      % record, by looking at volume history, does the volume increase ?
                % update iC, iE
                iC_current = struCell{iE_current}(iS).postCluster(iC_current);
                iE_current = iE_current + 1;
            else
                iC_current = 0;
                iE_current = iE_current + 1;
            end
            while iC_current
                struCell{iE_current}(iS).volEvo(iC_current,:) = vol;      % update volume in the cluster in the next iE
                struCell{iE_current}(iS).cvInc(iC_current) = vInc;
                % update iC, iE
                iC_current = struCell{iE_current}(iS).postCluster(iC_current);
                iE_current = iE_current + 1;
            end
            
            
        end
        

    end
end


%% (3) update to save
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save([saveDataPath,fName_c2t_result],'stru','-append');
end

