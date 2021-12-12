% chenzhe, 2019-07-17
% It is likely useful to manually label twin-gb intersection.  We need some ground truth data anyway.


clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the results where twin identification was based on trace dir and strain');

[twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the results for twin-grain boundary intersection');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
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
%% [data] strain data. Convert into v7.3 for partial loading
%%
clear strainFile;
for iE = iE_start-1:iE_stop
    strainFileName = [dicPath,'\',f2,STOP{iE+B}];
    disp(strainFileName);
    if ~exist([strainFileName,'_v73.mat'],'file')
        load(strainFileName);
        clear('exy_corrected');
        load(strainFileName,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
        
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
        
        outlier_removed = 1;
        save([strainFileName,'_v73.mat'],'outlier_removed','exy_corrected','-v7.3');
        
        myFile = matfile(strainFileName);
        myFields = who(myFile);
        for ii=1:length(myFields)
            save([strainFileName,'_v73.mat'],myFields{ii},'-append','-v7.3');
        end
    else
        disp('v7.3 file already exist');
    end
    strainFile{iE} = matfile([strainFileName,'_v73.mat']);
end
%% (0) load data, using SF threshold values to assign active twin system, and make maps
% Load cluster number maps (cleaned).

clusterNumberMapCell = cell(1,length(STOP)-1);
for iE = []%iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
% load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');
load(fullfile(confirmedLabelPath,confirmedLabelFile),'trueTwinMapCell');

% load previous twin_gb interaction result, for reference.
load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile));

%% Find triple points
[boundary,~, ~, tripleTF, ~, indTriple, tripleIDs] = find_one_boundary_from_ID_matrix(ID);
xTriple = X(indTriple);
yTriple = Y(indTriple);

%% Load struCell for updated tGbNormals here.

% Select an iE.
iE = 3;
iS = 1;
% get twinned grains from previous result
twinned_grain_list = [];
for iS = 1:length(struCell{iE})
    if sum(struCell{iE}(iS).cTrueTwin(:))>0
        twinned_grain_list = [twinned_grain_list; struCell{iE}(iS).gID];
    end
end

%% Select an iS to start
iS = 1;
variableNames = {'iE','ID','gDia','ID_neighbor','gDia_neighbor','TS','TSF','incoming',...
    'mPrime','rank_mPrime','ssn_neighbor','SF_neighbor',...
    'resB','rank_resB','ssn_neighbor_r','SF_neighbor_r',...
    'intersection_to_triple','twin_at_triple','first_appear_here','initiating','init_iE'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;
plotTF = 1;
pauseTF = 1;
T_template = T;

%% The purpose of this code is mainly compare intiation vs existing twin-slip pair, study the [effect of mPrime]
% Now I am thinking, do we need to consider intersection to triple, twin at triple, etc, for both sides of the grain boundry??????????????????????????????????   2019-08-27.   
warning('off','MATLAB:table:RowsAddedExistingVars');
continueTF = true;
while (continueTF)&&(iS<length(struCell{iE}))
    
    if sum(struCell{iE}(iS).cTrueTwin(:))>0  % ismember(struCell{iE}(iS).gID, twinned_grain_list)
        %%
        close all;
        ID_current = struCell{iE}(iS).gID;
        ind = find(gID==ID_current);
        
        euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        if (1==eulerAligned)
            % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
        else
            % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
            [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
        end
        [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
        traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
        
        % find out grain size
        gDia = sqrt(4*(struCell{iE}(iS).gVol* (360/4096*5)^2)/pi);
        
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);
        
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        X_local = X(indR_min:indR_max, indC_min:indC_max);
        Y_local = Y(indR_min:indR_max, indC_min:indC_max);
        uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
        boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        trueTwinMapLocal = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        e_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
        
        if plotTF==1
            %         [handleFig0,a1,~] = myplot(X_local, Y_local, e_local, boundaryTF_local);
            %         disableDefaultInteractivity(a1);
            %         [handleFig,aa,~] = myplot(X_local, Y_local, trueTwinMapLocal, boundaryTF_local);
            %         caxis([18 24]);
            [handleFig,aa,~] = myplot(X_local, Y_local, e_local, boundaryTF_local);  % boundaryTF_local_2 = find_boundary_from_ID_matrix(trueTwinMapLocal>0)|(boundaryTF_local);
            label_map_with_ID(X_local, Y_local, ID_local, handleFig, ID_current);
            disableDefaultInteractivity(aa);
            hold on;
        end
        
        % Find active system, if any, using cTrueTwin
        activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        tSF = struCell{iE}(iS).tSF;
        
        mPrime_table = [];
        resB_table = [];
        T_local = T_template;
        
        % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
        for iNb = 1:nNeighbors
            ID_neighbor = ID_neighbors(iNb);
            
            % (1.1) Calculate this_uniqueGB number.
            if ID_current > ID_neighbor
                this_uniqueGB = ID_current * 10000 + ID_neighbor;
            else
                this_uniqueGB = ID_neighbor * 10000 + ID_current;
            end
            
            % (1.2) Find which of the active twin systems in the grain are 'incoming' to this grain boundary
            % In addition, we need to find out if this twin first appear at this grain boundary.
            incoming_TS = zeros(1,6);
            ts_first_appear_here = zeros(1,6);  % if the twin was shown here at the iE that it just activated
            init_iE = zeros(1,6);   % the first iE that this twin shown at this gb  
            for iTwin = 1:6
                if ismember(this_uniqueGB, struCell{iE}(iS).tGb{iTwin})
                    incoming_TS(iTwin) = 1;
                    
                    % Find the first iE level that this twin was detected to intersect a grain boundary.
                    % If at that time, this boundary was intersected, then the ts first appear here.
                    ie=2;
                    while(ie<=iE)
                        if ~isempty(struCell{ie}(iS).tGb{iTwin})
                            break;
                        else
                            ie = ie + 1;
                        end
                    end
                    if ismember(this_uniqueGB, struCell{ie}(iS).tGb{iTwin})
                        ts_first_appear_here(iTwin) = 1;
                    end
                end
            end
            % incoming_TS_n = find(incoming_TS) + 18;
            
            if sum(incoming_TS)==0
                have_twin_intersect = 0;
            else
                have_twin_intersect = 1;
            end
            
            % [] check if incoming twin system intersect a triple point
            ind_of_indTriple = sum(ismember(tripleIDs,ID_current) + ismember(tripleIDs,ID_neighbor),2)==2;
            x_triple = X(indTriple(ind_of_indTriple));
            y_triple = Y(indTriple(ind_of_indTriple));
            
            twin_at_triple = zeros(6,1);
            intersection_to_triple = zeros(6,1);
            for iTwin = 1:6
                % If this twin variant intersects this grain boundary
                if incoming_TS(iTwin)==1
                    ind_tb = find(this_uniqueGB==struCell{iE}(iS).tGb{iTwin});      % find ind, because a twin can intersect several grain boundaries
                    tGbPts =struCell{iE}(iS).tGbPts{iTwin}{ind_tb};                 % also, this variant can intersect this grain boundary at several points
                    d_to_triple = pdist2(tGbPts,[x_triple,y_triple]);
                    [dmin,ind] = min(d_to_triple,[],2);     % all the points to its closest triple point
                    [dmax, ind_of_ind] = max(dmin,[],1);                  % the largest among distance_to_the_closest_triple_point
                    dmax = dmax * 360/4096;
                    intersection_to_triple(iTwin) = dmax;
                    % If this is still very small, say less than [3 micron] to triple point
                    if dmax < 3
                        twin_at_triple(iTwin) = 1;
                        text(x_triple(ind(ind_of_ind)),y_triple(ind(ind_of_ind)),50, [num2str(dmax)], 'color','b','fontsize',36);
                    end
                    
                    init_iE(iTwin) = min(struCell{iE}(iS).tGbPtsiE{iTwin}{ind_tb});
                    
                end
            end
            
            % ---------------> What if neighbor intersect close to a triple point?  Need to add some considerations?      ---------------------------???????????    
            % (1.3) Similarly, find if this neighbor grain (ID_neighbor) is twinned.
            % If twinned, and twin intersect this_uniqueGB, consider only these intersecting twins
            % else, consider only basal in this_neighbor
            iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
            if ~isempty(iS_neighbor)
                gDia_neighbor = sqrt(4*(struCell{iE}(iS_neighbor).gVol* (360/4096*5)^2)/pi);
                activeTS_Nb = sum(struCell{iE}(iS_neighbor).cTrueTwin,1)>0;
                for iTwin = 1:6
                    if ~ismember(this_uniqueGB, struCell{iE}(iS_neighbor).tGb{iTwin})
                        activeTS_Nb(iTwin) = 0;
                    end
                end
                if sum(activeTS_Nb)==0
                    % Only basal in this_Nb
                    outgoing_basal = 1;
                    outgoing_ss_n = [1 2 3];
                else
                    outgoing_basal = 0;
                    outgoing_ss_n = 19:24;
                    outgoing_ss_n = outgoing_ss_n(logical(activeTS_Nb));
                end
                
                % Calculate mPrime
                ind_Nb = find(gID==ID_neighbor);
                euler_Nb = [gPhi1(ind_Nb),gPhi(ind_Nb),gPhi2(ind_Nb)];
                [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler, euler_Nb, stressTensor, [1 0 0], sampleMaterial, 'twin');
                
                [max_basal_SF_Nb, ~] = max(schmidFactorG2(1:3));
                [max_twin_SF_Nb, ~] = max(schmidFactorG2(19:24));
                
                % Based on whether considering basal or twin in this neighbor, do different things: ...
                % (3) label basal/twin slip trace
                ID_target = ID_neighbor;
                if plotTF==1
                    local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, outgoing_ss_n, gca);
                end
                % find max basal SF, and the corresponding basal system, in the neighbor grain
                SFs_Nb = schmidFactorG2(outgoing_ss_n);
                
                
                % [[[[--> maybe it's better to find potential basal. e.g., all those with sf > 0.3, or normalized > 0.8
                % [max_sf, ind_max_sf] = max(SFs_Nb);
                potential_ss = (SFs_Nb>0.25) | (SFs_Nb./max(SFs_Nb)>0.8);   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                
                mPrime_local = mPrimeMatrixAbs(19:24, outgoing_ss_n);
                [mPrime_each_twin, ind_in_outgoing_n] = max(mPrime_local.*reshape(potential_ss,1,[]),[],2);
                ssn_each_twin = outgoing_ss_n(ind_in_outgoing_n);
                SF_each_twin = SFs_Nb(ind_in_outgoing_n);
                mPrime_table = [mPrime_table, mPrime_each_twin];
                [sorted_mPrime_each_twin,~] = sort(mPrime_each_twin,'descend');
                [~,rank_mPrime_each_twin] = ismember(mPrime_each_twin, sorted_mPrime_each_twin);
                % find the mPrime to plot
                if have_twin_intersect
                    [mPrime, ind_mPrime] = max(mPrime_each_twin.*incoming_TS(:));
                    mPrime_rank = rank_mPrime_each_twin(ind_mPrime);
                else
                    [mPrime, ind_mPrime] = max(mPrime_each_twin.*activeTS(:));
                    mPrime_rank = rank_mPrime_each_twin(ind_mPrime);
                end
                
                
                resB_local = resBurgersMatrixAbs(19:24, outgoing_ss_n);
                tt = resB_local.*reshape(potential_ss,1,[]);
                tt = tt./logical(tt);
                [resB_each_twin, ind_in_outgoing_n] = min(tt,[],2);
                ssn_each_twin_r = outgoing_ss_n(ind_in_outgoing_n);
                SF_each_twin_r = SFs_Nb(ind_in_outgoing_n);
                resB_table = [resB_table, resB_each_twin];
                [sorted_resB_each_twin,~] = sort(resB_each_twin,'ascend');
                [~,rank_resB_each_twin] = ismember(resB_each_twin, sorted_resB_each_twin);
                % find the resB to plot
                if have_twin_intersect
                    [resB, ind_resB] = min(resB_each_twin./incoming_TS(:));
                    resB_rank = rank_resB_each_twin(ind_resB);
                else
                    [resB, ind_resB] = min(resB_each_twin./activeTS(:));
                    resB_rank = rank_resB_each_twin(ind_resB);
                end
                
                % [NOTE] --> Draw mPrime and color, but label ResB, to maximize information provided.
                if plotTF==1
                    inds = find(uniqueBoundary_local==this_uniqueGB);
                    plot(X_local(inds),Y_local(inds),'.','markersize',18,'color',[1 1 0] * mPrime);
                    if outgoing_basal==1
                        text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_rank + mPrime,4)], 'color','r','fontsize',16);
                    else
                        text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_rank + mPrime,4)], 'color','b','fontsize',16);
                    end
                end
                
                % Extend table for summary.  -- can assign in local table using name?
                for iTwin=1:6
                    % T_local = [T_local;{iE, ID_current, gDia, ID_neighbor, gDia_neighbor, iTwin, tSF(iTwin), incoming_TS(iTwin), ...
                    %     mPrime_each_twin(iTwin), rank_mPrime_each_twin(iTwin), ssn_each_twin(iTwin), SF_each_twin(iTwin), ...
                    %     resB_each_twin(iTwin), rank_resB_each_twin(iTwin), ssn_each_twin_r(iTwin), SF_each_twin_r(iTwin), ...
                    %     intersection_to_triple(iTwin), twin_at_triple(iTwin), ts_first_appear_here(iTwin), false}];
                    ir = size(T_local,1) + 1;
                    T_local.iE(ir) = iE;
                    T_local.ID(ir) = ID_current;
                    T_local.gDia(ir) = gDia;
                    T_local.ID_neighbor(ir) = ID_neighbor;
                    T_local.gDia_neighbor(ir) = gDia_neighbor;
                    T_local.TS(ir) = iTwin;
                    T_local.TSF(ir) = tSF(iTwin);
                    T_local.incoming(ir) = incoming_TS(iTwin);
                    T_local.mPrime(ir) = mPrime_each_twin(iTwin);
                    T_local.rank_mPrime(ir) = rank_mPrime_each_twin(iTwin);
                    T_local.ssn_neighbor(ir) = ssn_each_twin(iTwin);
                    T_local.SF_neighbor(ir) = SF_each_twin(iTwin);
                    T_local.resB(ir) = resB_each_twin(iTwin);
                    T_local.rank_resB(ir) = rank_resB_each_twin(iTwin);
                    T_local.ssn_neighbor_r(ir) = ssn_each_twin_r(iTwin);
                    T_local.SF_neighbor_r(ir) = SF_each_twin_r(iTwin);
                    T_local.intersection_to_triple(ir) = intersection_to_triple(iTwin);
                    T_local.twin_at_triple(ir) = twin_at_triple(iTwin);
                    T_local.first_appear_here(ir) = ts_first_appear_here(iTwin);
                    T_local.initiating(ir) = false;
                    T_local.init_iE(ir) = init_iE(iTwin);
                end
            end
        end
        
        % 'incoming_TS' means intersecting.  Each twin can intersect several gbs. Find the one with the (highest m') to represent the initiating point.
        for iTwin = 1:6
            rows = find((T_local.TS==iTwin)&(T_local.first_appear_here==1));
            if ~isempty(rows)
                [~,ind] = max(T_local.mPrime(rows,:));
                T_local.initiating(rows(ind)) = 1;
                disp(['Twin ',num2str(iTwin),' first appear, initiatiated from neighbor ID:',num2str(T_local.ID_neighbor(rows(ind)))]);
            end
        end
        ind = (T_local.initiating==1);
        label_map_with_ID(X_local,Y_local,ID_local,gcf,T_local.ID_neighbor(ind),'m');
        T = [T;T_local];
        %%
        if plotTF==1
            % (3) Draw previously labeled lines
            for iTwin = 1:6
                tGbNormals = struCell{iE}(iS).tGbNormal{iTwin};
                for iGb = 1:length(tGbNormals)
                    if ~isempty(tGbNormals{iGb})
                        for iLine = 1:length(tGbNormals{iGb})
                            drawline(aa,'Position',tGbNormals{iGb}{iLine},'color','m','linewidth',1);
                        end
                    end
                end
            end
            
            if pauseTF==1
                needToAdd = questdlg('Continue, pause, or cancel','select answer','continue','pause','cancel','continue');
                switch needToAdd
                    case 'pause'
                        continueTF = false;
                        hgexport(gcf,'-clipboard');
                    case 'continue'
                        continueTF = true;
                        hgexport(gcf,'-clipboard');
                    case 'cancel'
                        iS = iS - 1;    % when cancel, reduce by 1 first
                        continueTF = false;
                end
            end
            saveas(gcf,['temp_results\iE_',num2str(iE),'_ID_',num2str(ID_current),'.tiff']);
        end
    end
    disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
    iS = iS + 1;
end

warning on;
%% it's better to save periodically
try
    save(['temp_results\for_recover_iE_',num2str(iE)],'struCell','iE','iS','T','-append');
catch
    save(['temp_results\for_recover_iE_',num2str(iE)],'struCell','iE','iS','T');
end

%% Record data at all iEs
iE_current = iE;    % record

% initiate
tb_gbNum = [];
tb_iE = [];
tb_gNum = [];
tb_tsNum = [];
tb_pts = [];
tb_iE_inGrain = []; % first iE level appear in this grain
tb_iE_atBoundary = [];  % first iE level appear at this boundary

tboundary = zeros(size(uniqueBoundary));

for iE = iE_start:iE_current
    % summarize
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID
        activeTS = find(sum(struCell{iE}(iS).cTrueTwin, 1)> 0);     % but maybe haven't intersected a boundary yet
        for ii = 1:length(activeTS)
            iTwin = activeTS(ii);
            tsNum = iTwin + 18;
            
            % find unique grain id, grain ID, etc
            GBs = struCell{iE}(iS).tGb{iTwin};
            Pts = struCell{iE}(iS).tGbPts{iTwin};
            iEs = struCell{iE}(iS).tGbPtsiE{iTwin};
            
            iE_inGrain = min(cell2mat(struCell{iE}(iS).tGbPtsiE{iTwin}(:)));
            
            for iGb = 1:length(GBs)
                iE_atBoundary = min(struCell{iE}(iS).tGbPtsiE{iTwin}{iGb});
                
                ind = find((tb_iE==iE)&(tb_gbNum==GBs(iGb))&(tb_gNum==ID_current)&(tb_tsNum==tsNum));
                if isempty(ind)
                    ind = size(tb_gbNum,1)+1;
                    tb_iE(ind,1) = iE;
                    tb_gbNum(ind,1) = GBs(iGb);
                    tb_gNum(ind,1) = ID_current;
                    tb_tsNum(ind,1) = tsNum;
                    tb_pts{ind,1} = Pts{iGb};
                    
                    tb_iE_inGrain(ind,1) = iE_inGrain;
                    tb_iE_atBoundary(ind,1) = iE_atBoundary;
                    tboundary(uniqueBoundary == GBs(iGb)) = 1;
                else
                    disp('error');
                end
            end
        end
    end
    
    vN = {'iE','gb','g','ts','iE_g','iE_gb'};
    % TA = cell2table(cell(0,length(vN)));
    TA = array2table([tb_iE, tb_gbNum, tb_gNum, tb_tsNum, tb_iE_inGrain, tb_iE_atBoundary]);
    TA.Properties.VariableNames = vN;
    
    % plot, be careful
    tBoundaryCell{iE} = [];
    tBoundaryCell{iE} = tboundary;
end

iE = iE_current;    % change back

%%
timeStr = datestr(now,'yyyymmdd_HHMM');
save(['temp_results\',timeStr,'_iE_',num2str(iE),'_twin_at_boundary_manual_result.mat'], 'TA', 'tb_gbNum', 'tb_iE', 'tb_gNum', 'tb_tsNum', 'tb_pts', 'tb_iE_inGrain', 'tb_iE_atBoundary',...
    'tBoundaryCell','tBoundary_accum','struCell','-v7.3');
% save([timeStr,'_twin_at_boundary_result_ws.mat'],'-v7.3');

%% Summary. To determine twin activation, effect of (1) m' rank, (2) m', (3) resB, (4) resB rank
% Objects to compare: twinned variants in twinned grains
% (1) All twin boundaries
% (2) Boundaries with intersecting twins
% (3) Boundaries where twins are considered to be initiated from
% Properties to compare:
% (1) m', (2) m' rank, (3) resB, (4) resB rank

%% (1) m'
close all;
% (1) all boundaries, distribution of m'
figure;
histogram(T.mPrime, 0:0.1:1);
xlabel('m''');
ylabel('Counts');
title('all boundaries, all twins');

% (2) boundaries with intersecting twins
ind = (T.incoming==1);
t = T(ind,:);
figure;
histogram(t.mPrime, 0:0.1:1);
xlabel('m''');
ylabel('Counts');
title('only intersecting twins');

% (3) boundaries with initiating twins
ind = (T.incoming==1)&(T.initiating==1);
t = T(ind,:);
figure;
histogram(t.mPrime, 0:0.1:1);
xlabel('m''');
ylabel('Counts');
title('initiating twins');
ylim = get(gca,'ylim');

% (4) boundaries with initiating twins that are not at triple points
ind = (T.incoming==1)&(T.initiating==1)&(T.twin_at_triple==0);
t = T(ind,:);
figure;
histogram(t.mPrime, 0:0.1:1);
xlabel('m''');
ylabel('Counts');
title('initiating twins not at triple points');
set(gca,'ylim',ylim);

%% (2) m' rank
close all;
% (1) all boundaries, distribution of m'
figure;
histogram(T.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('all boundaries, all twins');

% (2) boundaries with intersecting twins
ind = (T.incoming==1);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('only intersecting twins');

% (3) boundaries with initiating twins
ind = (T.incoming==1)&(T.initiating==1);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('initiating twins');
ylim = get(gca,'ylim');

% (4) boundaries with initiating twins that are not at triple points
ind = (T.incoming==1)&(T.initiating==1)&(T.twin_at_triple==0);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('initiating twins not at triple points');
set(gca,'ylim',ylim);

% % (5) divide by grain size
% ind = (T.incoming==1)&(T.initiating==1)&(T.twin_at_triple==0)&(T.gDia<150);
% t = T(ind,:);
% figure;
% histogram(t.rank_mPrime, 0.5:6.5);
% xlabel('m'' rank');
% ylabel('Counts');
% title('initiating twins not at triple points');
% set(gca,'ylim',ylim);


%% (3) resB
close all;
% (1) all boundaries, distribution of m'
figure;
histogram(T.resB, 0:0.1:1);
xlabel('resB');
ylabel('Counts');
title('all boundaries, all twins');

% (2) boundaries with intersecting twins
ind = (T.incoming==1);
t = T(ind,:);
figure;
histogram(t.resB, 0:0.1:1);
xlabel('resB');
ylabel('Counts');
title('only intersecting twins');

% (3) boundaries with initiating twins
ind = (T.incoming==1)&(T.initiating==1);
t = T(ind,:);
figure;
histogram(t.resB, 0:0.1:1);
xlabel('resB');
ylabel('Counts');
title('initiating twins');
ylim = get(gca,'ylim');

% (4) boundaries with initiating twins that are not at triple points
ind = (T.incoming==1)&(T.initiating==1)&(T.twin_at_triple==0);
t = T(ind,:);
figure;
histogram(t.resB, 0:0.1:1);
xlabel('resB');
ylabel('Counts');
title('initiating twins not at triple points');
set(gca,'ylim',ylim);

%% (4) resB
close all;
% (1) all boundaries, distribution of m'
figure;
histogram(T.rank_resB, 0.5:6.5);
xlabel('resB rank');
ylabel('Counts');
title('all boundaries, all twins');

% (2) boundaries with intersecting twins
ind = (T.incoming==1);
t = T(ind,:);
figure;
histogram(t.rank_resB, 0.5:6.5);
xlabel('resB rank');
ylabel('Counts');
title('only intersecting twins');

% (3) boundaries with initiating twins
ind = (T.incoming==1)&(T.initiating==1);
t = T(ind,:);
figure;
histogram(t.rank_resB, 0.5:6.5);
xlabel('resB rank');
ylabel('Counts');
title('initiating twins');
ylim = get(gca,'ylim');

% (4) boundaries with initiating twins that are not at triple points
ind = (T.incoming==1)&(T.initiating==1)&(T.twin_at_triple==0);
t = T(ind,:);
figure;
histogram(t.rank_resB, 0.5:6.5);
xlabel('resB rank');
ylabel('Counts');
title('initiating twins not at triple points');
set(gca,'ylim',ylim);

%% (5) m' vs neighbor_SF

