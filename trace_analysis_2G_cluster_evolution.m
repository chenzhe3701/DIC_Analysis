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

%% select iEs to analyze cluster evolution
cNumMaps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0 
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru');    
    cNumMaps{iE} = clusterNumMap;
    struCell{iE} = stru;
end

%%
for iE = iE_start:iE_stop-1
    %%
    %     iE = 4;
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    hWaitbar = waitbar(0,'running each grain, each cluster, ...');
    for iS =1:length(struA)
        % iS = find(arrayfun(@(x) x.gID == 246,stru));  % for debugging
        % iS = find(gIDwithTrace == 296); % for debugging.
        % iS = find(arrayfun(@(x) x.gID == 190,stru));
        
        ID_current = gIDwithTrace(iS);

        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        
        cMapA = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
        cMapP = cNumMaps{iE+1}(indR_min:indR_max, indC_min:indC_max);
        
        cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
        cMapP(ID_local~=ID_current) = 0;  % cluster number just this grain
        
        if ~isfield(struA(iS),'vrFwd')||isempty(struA(iS).vrFwd)
            struA(iS).vrFwd = zeros(size(struA(iS).cLabel));
        end
        if ~isfield(struA(iS),'vrBwd')||isempty(struA(iS).vrBwd)
            struA(iS).vrBwd = zeros(size(struA(iS).cLabel));
        end
        if ~isfield(struP(iS),'vrFwd')||isempty(struP(iS).vrFwd)
            struP(iS).vrFwd = zeros(size(struP(iS).cLabel));
        end
        if ~isfield(struP(iS),'vrBwd')||isempty(struP(iS).vrBwd)
            struP(iS).vrBwd = zeros(size(struP(iS).cLabel));
        end
        
        % match clusters, and find the one with area grown.
        cOverlap = [];
        for ii=1:length(struA(iS).cLabel)
            cNumA = struA(iS).cLabel(ii);
            for jj=1:length(struP(iS).cLabel)
                cNumP = struP(iS).cLabel(jj);
                cOverlap(ii,jj) =sum(sum((cMapA==cNumA)&(cMapP==cNumP)));
            end
        end
        volA = sum(cOverlap,2);
        volP = sum(cOverlap,1);
        vrFwd = volP./volA;
        vrBwd = volA./volP;
        overlapPctA = cOverlap./volA;
        [cFrom,cTo] = hungarian_assign(max(cOverlap(:))-cOverlap);
        link = false(size(cOverlap));
        for ii = 1:length(cFrom)
            if (cFrom(ii)>0)&&(cTo(ii)>0)
                link(cFrom(ii),cTo(ii)) = true;
            end
        end
        twinLikely = link & (vrFwd<=1) & (overlapPctA>0.5);
        [inda,indp] = find(twinLikely);
        struA(iS).vrFwd(inda) = vrFwd(indp);
        struP(iS).vrBwd(indp) = vrBwd(inda);
        
%         myplot(cMapA);
%         myplot(cMapP);
        waitbar(iS/length(stru), hWaitbar);
    end
    
    try
        close(hWaitbar);
    end
    warning('on','all');
    
    struCell{iE} = struA;   % pre
    struCell{iE+1} = struP; % post
    
end

%% update to save

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    stru = struCell{iE};
    save([saveDataPath,fName_c2t_result],'stru','-append');    
    cNumMaps{iE} = clusterNumMap;
    struCell{iE} = stru;
end


