% track cluster NumEvolution and VolEvolution
% chenzhe, 2018-03-01


% clear;
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

modify_cVol_TF = 0; % if already have 'cVol', do not modify again. very important!
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
        stru = rmfield(stru,{'tProbAvg'});
    end
    % initialize/zero fields
    for iS =1:length(stru)
        if modify_cVol_TF
            stru(iS).cVol = zeros(size(stru(iS).cLabel));   % the volume of this cluster
            stru(iS).cVolCleaned = zeros(size(stru(iS).cLabel)); % the cluster volume overlap with postCluster, and cleaned
        end
        stru(iS).preCluster = zeros(size(stru(iS).cLabel));
        stru(iS).postCluster = zeros(size(stru(iS).cLabel));
        
        stru(iS).volEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).volEvoCleaned = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbMax = zeros(size(stru(iS).cLabel));
        stru(iS).cvInc = zeros(size(stru(iS).cLabel));      % based on looking at the volume evolution in all strain levels
    end
    struCell{iE} = stru;
end

% (1) match clusters and track volume evolution, fields related: cVol, cVolCleaned, preCluster, postCluster
for iE = iE_start:iE_stop-1
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(struA)
        %% iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        % iS = find(gIDwithTrace == 296); % for debugging.
        %kk=kk+1,close,close
        %iS = find(arrayfun(@(x) x.gID == varEnabled(kk,1),stru));
        
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
                
                if modify_cVol_TF
                    struA(iS).cVol(ii) = sum(cMapA(:)==cNumA);  % count the volume for the current cluster
                    struA(iS).cVolCleaned(ii) = sum(cMapA_cleaned(:)==cNumA);  % count the volume for the current cluster
                    if iE==(iE_stop-1)
                        struP(iS).cVol(jj) = sum(cMapP(:)==cNumP);  % count the volume for the postCluster
                        struP(iS).cVolCleaned(jj) = sum(cMapP(:)==cNumP);  % count the volume for the postCluster
                    end
                end
                cOverlap(ii,jj) =sum(sum((cMapA_cleaned==cNumA)&(cMapP_cleaned==cNumP)));   % use cleaned map to calculate how good they overlap !!! But NOT the actual cluster size ------------------
            end
        end
        volA_cleaned = sum(cOverlap,2);
            
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
        
        %         myplot(cMapA_cleaned);
        %         myplot(cMapP_cleaned);
        %         run('D:\p\m\DIC_Analysis\tempfun\script_help_summarize_2018_03_04.m')
        
        
    end
    
    try
        close(hWaitbar);
    end
    warning('on','all');
    
    struCell{iE} = struA;   % pre
    struCell{iE+1} = struP; % post
    
end

% (after 1) update to Save. Can disable to temporarily to prevent data loss.
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save([saveDataPath,fName_c2t_result],'stru','-append');
end


%% (2) after tracking all strain level pairs, get volume history in the whole strain levels, and analyze. Related fields: volEvo, volEvoCleaned, cvInc
%  First, re-load stru in all strain levels
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    try
        stru = rmfield(stru,{'tProbAvg'});
    end
    % initialize/zero related fields
    for iS =1:length(stru)
        stru(iS).volEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).volEvoCleaned = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbEvo = zeros(length(stru(iS).cLabel),length(STOP)-1);
        stru(iS).tProbMax = zeros(size(stru(iS).cLabel));
        stru(iS).cvInc = zeros(size(stru(iS).cLabel));      % based on looking at the volume evolution in all strain levels
    end
    struCell{iE} = stru;
end


for iE = iE_start:iE_stop
    % stru = struCell{iE};
    for iS =1:length(struCell{iE})
        
        for iCluster = 1:length(struCell{iE}(iS).cLabel)
            
            vol = zeros(1,length(STOP)-1);
            vol_cleaned = zeros(1,length(STOP)-1);
            
            % [1st loop] find the volume of this area in different strain levels
            iE_list = iE;
            iC_list = iCluster;
            vol = struCell{iE}(iS).cVol(iCluster);    % record size of the current cluster in the current iE
            vol_cleaned = struCell{iE}(iS).cVolCleaned(iCluster);    % record size of the current cluster overlaid with the post cluster, and cleaned
            tProbEvo = struCell{iE}(iS).tProb(iCluster);
            
            % search to earlier strain (or try not to...)
            while 0 ~= struCell{iE_list(1)}(iS).preCluster(iC_list(1))
                iC_list = [struCell{iE_list(1)}(iS).preCluster(iC_list(1)), iC_list];
                iE_list = [iE_list(1)-1, iE_list];
                vol = [struCell{iE_list(1)}(iS).cVol(iC_list(1)), vol];
                vol_cleaned = [struCell{iE_list(1)}(iS).cVolCleaned(iC_list(1)), vol_cleaned];
                tProbEvo = [struCell{iE_list(1)}(iS).tProb(iC_list(1)), tProbEvo];
            end
            % search to later strain
            while 0 ~= struCell{iE_list(end)}(iS).postCluster(iC_list(end))
                iC_list = [iC_list,struCell{iE_list(end)}(iS).postCluster(iC_list(end))];
                iE_list = [iE_list, iE_list(end)+1];
                vol = [vol, struCell{iE_list(end)}(iS).cVol(iC_list(end))];
                vol_cleaned = [vol_cleaned, struCell{iE_list(end)}(iS).cVolCleaned(iC_list(end))];
                tProbEvo = [tProbEvo, struCell{iE_list(end)}(iS).tProb(iC_list(end))];
            end
            
            
            % evaluate if volume seems to be increasing. Is this method OK?
            vol_valid = vol_cleaned;   % volume evolution in whole history. Use 'vol' or use 'vol_cleaned'?
            if 1==length(vol_valid)
                cvInc = 0;
            else
                cvInc = diff(vol_valid)./conv(vol_valid, [0.5, 0.5], 'valid');
                cvInc = mean(cvInc);
            end
            
            % [2nd loop] copy
            if length(iE_list) >= 1
                for ii = 1:length(iE_list)
                    iE_to_assign = iE_list(ii);
                    iC_to_assign = iC_list(ii);
                    
                    struCell{iE_to_assign}(iS).volEvo(iC_to_assign,iE_list) = vol;    % record size history of the current cluster in the current iE
                    struCell{iE_to_assign}(iS).volEvoCleaned(iC_to_assign,iE_list) = vol_cleaned;
                    struCell{iE_to_assign}(iS).tProbEvo(iC_to_assign,iE_list) = tProbEvo;
                    
                    struCell{iE_to_assign}(iS).cvInc(iC_to_assign) = cvInc;
                    
                    % The search is form iE_start to iE_to_assign
                    struCell{iE_to_assign}(iS).tProbMax(iC_to_assign) = max(tProbEvo(1:find(iE_list==iE_to_assign)));
                end
                
            end
            
            
        end
        
        
    end
end

% (after 2) update to save. Can disable to temporarily to prevent data loss.
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save([saveDataPath,fName_c2t_result],'stru','-append');
end


%% [for debug directly]
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    cNumMaps{iE} = clusterNumMap;
    cNumMaps_cleaned{iE} = clusterNumMapCleaned;
    
    struCell{iE} = stru;
end
%% [for debug] plot map of a cluster of intereted in all strain levels

egc = 211442;
egc = 210121;
iC_target = mod(egc,10);
egc = (egc-iC_target)/10;
ID_target = mod(egc,10000);
iE_target = (egc-ID_target)/10000;


debug = 1;
if debug
    close all;

    iS = find(arrayfun(@(x) x.gID == ID_target,struCell{iE_start}));
    
    % a good method to, from the starting iE and iC, fill all useful iEs and iCs
    iE_list = iE_target;
    iC_list = iC_target;
    while 0 ~= struCell{iE_list(1)}(iS).preCluster(iC_list(1))
        iC_list = [struCell{iE_list(1)}(iS).preCluster(iC_list(1)), iC_list];
        iE_list = [iE_list(1)-1, iE_list];
    end
    
    while 0 ~= struCell{iE_list(end)}(iS).postCluster(iC_list(end))
        iC_list = [iC_list,struCell{iE_list(end)}(iS).postCluster(iC_list(end))];
        iE_list = [iE_list, iE_list(end)+1];
    end
    
    
    ID_current=ID_target;
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    
    for ii = 1:length(iE_list)
        
        iE = iE_list(ii);
        iC = iC_list(ii);
        
        cMapA = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapA_cleaned = cNumMaps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
        
        cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
        cMapA_cleaned(ID_local~=ID_current) = 0;
        
        cMapA(cMapA==iC) = cMapA(cMapA==iC);
        myplot(cMapA);
        myplot(cMapA_cleaned);
        
    end
end

%%
close all
for iE = 2
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    
    iS = find(arrayfun(@(x) x.gID == ID_target,stru));  % for debugging
    
    
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
    cOverlap_noClean = [];
    cOverlap = [];
    for ii=1:length(struA(iS).cLabel)
        cNumA = struA(iS).cLabel(ii);
        for jj=1:length(struP(iS).cLabel)
            cNumP = struP(iS).cLabel(jj);
            %%%%%%%%%%%%%%%
            cOverlap(ii,jj) =sum(sum((cMapA_cleaned==cNumA)&(cMapP_cleaned==cNumP)));   % use cleaned map to calculate how good they overlap !!! But NOT the actual cluster size ------------------
            cOverlap_noClean(ii,jj) =sum(sum((cMapA==cNumA)&(cMapP==cNumP)));
        end
    end
    volA_cleaned = sum(cOverlap,2);
    
    % this is more precise:
    volA_noClean = struA(iS).cVol;
    volA_cleaned = struA(iS).cVolCleaned;
    
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
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    myplot(cMapA_cleaned);
    myplot(cMapP_cleaned);

    struCell{iE} = struA;   % pre
    struCell{iE+1} = struP; % post
    
end
run('D:\p\m\DIC_Analysis\tempfun\script_help_summarize_2018_03_04.m')