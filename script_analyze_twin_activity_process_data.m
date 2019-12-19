% chenzhe, 2019-09-03
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
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
% [confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the results where twin identification was based on trace dir and strain');

% [twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\*.mat','select the results for twin-grain boundary intersection');

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

% [data] strain data. Convert into v7.3 for partial loading
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

% (0) load data, using SF threshold values to assign active twin system, and make maps
% Load cluster number maps (cleaned).

clusterNumberMapCell = cell(1,length(STOP)-1);
for iE = []%iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
% load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');
% load(fullfile(confirmedLabelPath,confirmedLabelFile),'trueTwinMapCell');

% load previous twin_gb interaction result, for reference.
% load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile));

[newVariantFile, newVariantFilePath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the new result of dividing twin into variants');
load(fullfile(newVariantFilePath,newVariantFile),'struCell');   %,'trueTwinMapCell');

%% Find triple points
%%
[boundary,~, ~, tripleTF, ~, indTriple, tripleIDs] = find_one_boundary_from_ID_matrix(ID);
xTriple = X(indTriple);
yTriple = Y(indTriple);
%% Select an iS to start
%%
plotTF = 1;
pauseTF = 0;
umPerX = 360/4096;
umPerDp = 360/4096*5;

for iE = 2:5
    variableNames = {'iE','ID','gDia','ID_neighbor','gDia_neighbor','TS','TSF','gb_length','gb_dir',...
        'incoming','iiE_each_twin','iiE_each_twin_at_this_boundary','intersection_to_triple','iiE_twins_at_this_boundary_nb',...
        'mPrime','rank_mPrime','ssn_nb','SF_nb',...
        'resB','rank_resB','ssn_nb_r','SF_nb_r',...
        'mPrime_wrtB','rank_mPrime_wrtB','ssn_nb_wrtB','SF_nb_wrtB', 'resB_wrtB','rank_resB_wrtB','ssn_nb_r_wrtB','SF_nb_r_wrtB',...
        'mPrime_wrtT','rank_mPrime_wrtT','ssn_nb_wrtT','SF_nb_wrtT', 'resB_wrtT','rank_resB_wrtT','ssn_nb_r_wrtT','SF_nb_r_wrtT',...
        'initiating', 'eMean_1','eMean_2','eMean_1_nb','eMean_2_nb', 'max_basal_SF','max_twin_SF','max_basal_SF_nb','max_twin_SF_nb',...
        'exz_ba','exz_pr','exz_py','exz_pyII','exz_etw', 'exzr_ba','exzr_pr','exzr_py','exzr_pyII','exzr_etw',...
        'tGbVol','tGbVolPct','tGbStrength'};
    T = cell2table(cell(0,length(variableNames)));
    T.Properties.VariableNames = variableNames;
    T2 = T;
    T_template = T;
    lookupMa = [];
    lookupMb = [];
    
    eMap = calculate_effective_strain(strainFile{iE-1}.exx, strainFile{iE-1}.exy, strainFile{iE-1}.eyy);
    edmat = [];
    
    % summary of [bNum, gNum] pairs list, for different categories of activities
    % (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,
    % (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
    % (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth
    bg_not_involved = [];
    bg_slip_twin_a = [];
    bg_slip_twin_b =[];
    bg_co_found = [];
    bg_twin_twin_a = [];
    bg_twin_twin_b = [];
    bg_slip_growth_a = [];
    bg_slip_growth_b = [];
    bg_co_growth = [];
    
    % The purpose of this code is mainly compare intiation vs existing twin-slip pair, study the [effect of mPrime]
    iS = 1;
    warning('off','MATLAB:table:RowsAddedExistingVars');
    continueTF = true;
    dToTriple_th = 5;       % eliminate intersection whose distance to triple point is smaller than this value
    dToTriple_th_to_label = dToTriple_th;  % label if distance of intersection to triple point is smaller than this value
    
    hW = waitbar(0, ['iE=',num2str(iE),' analyze each grain']);
    hN = length(struCell{iE});
    while (continueTF)&&(iS<=length(struCell{iE}))
        waitbar(iS/hN, hW);
        
        close all;
        ID_current = struCell{iE}(iS).gID
        ind = find(gID==ID_current);
        
        euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
%         if (1==eulerAligned)
%             % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
%             [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
%         else
%             % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
%             [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
%         end
%         [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
%         traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
        
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
        exx_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);  % This is for effective strain
        
        % find out grain size, in um
        gDia = sqrt(4*(struCell{iE}(iS).gVol* (umPerDp)^2)/pi);
        % Find active system, if any, using cTrueTwin/tGb field
        % activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        activeTS = cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb);
        tSF = struCell{iE}(iS).tSF;
        
        % If grain of interest is twinned (here it means: got a twin-gb intersection labeled)
        if any(cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb)) % sum(struCell{iE}(iS).cTrueTwin(:))>0
            if plotTF==1
                %         [handleFig0,a1,~] = myplot(X_local, Y_local, e_local, boundaryTF_local);
                %         disableDefaultInteractivity(a1);
                %         [handleFig,aa,~] = myplot(X_local, Y_local, trueTwinMapLocal, boundaryTF_local);
                %         caxis([18 24]);
                [handleFig,aa,~] = myplot_bg(X_local, Y_local, exx_local, boundaryTF_local);  % boundaryTF_local_2 = find_boundary_from_ID_matrix(trueTwinMapLocal>0)|(boundaryTF_local);
                label_map_with_ID(X_local, Y_local, ID_local, handleFig, ID_current);
                disableDefaultInteractivity(aa);
                hold on;
            end
            
            T_local = T_template;
            
            % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
            for iNb = 1:nNeighbors
                ID_neighbor = ID_neighbors(iNb);
                iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
                
                ind = find(gID==ID_neighbor);
                euler_nb = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
                
                if ~isempty(iS_neighbor)
                    gDia_neighbor = sqrt(4*(struCell{iE}(iS_neighbor).gVol* (360/4096*5)^2)/pi);
                    % activeTS_nb = sum(struCell{iE}(iS_neighbor).cTrueTwin,1)>0;
                    activeTS_nb = cellfun(@(x) ~isempty(x), struCell{iE}(iS_neighbor).tGb);
                    
                    
                    % (1.1) Calculate this_uniqueGB number.
                    if ID_current > ID_neighbor
                        gb = ID_current * 10000 + ID_neighbor;
                    else
                        gb = ID_neighbor * 10000 + ID_current;
                    end
                    % Find gb_length and gb_dir
                    ind = uniqueBoundary_local(:)==gb;
                    gb_length = sum(ind);
                    xs = X_local(ind);
                    ys = Y_local(ind);
                    mdl = fitlm(xs,ys);
                    gb_dir = atand(mdl.Coefficients.Estimate(2));
                    gb_normal = 90 + gb_dir;
                    gb_normal_v = [cosd(gb_normal),sind(gb_normal),0]; 
                    gb_normal_v = gb_normal_v/norm(gb_normal_v);
                    
                    % strain calculation in area of interest.
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                    ind = (distMap_local>0)&(distMap_local<=50);
                    eMean_1 =  nanmean(eMap_local(ind));   % mean
                    ind = (distMap_local>50)&(distMap_local<=100);
                    eMean_2 =  nanmean(eMap_local(ind));   % mean
                    
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_neighbor]);
                    ind = (distMap_local>0)&(distMap_local<=50);
                    eMean_1_nb =  nanmean(eMap_local(ind));   % mean
                    ind = (distMap_local>50)&(distMap_local<=100);
                    eMean_2_nb =  nanmean(eMap_local(ind));   % mean
                    
                    ind_of_indTriple = sum(ismember(tripleIDs,ID_current)+ismember(tripleIDs,ID_neighbor), 2)==2;
                    xTriple_local = X(indTriple(ind_of_indTriple));
                    yTriple_local = Y(indTriple(ind_of_indTriple));
                    
                    % Find [iE, iTwin, dToTriple] on the grain_of_interest side, then on the neighboring side
                    [iE_iTwin_dT_list, valid_grain_a] = find_activity_at_boundary_from_struCell(struCell, iE, ID_current, gb, [xTriple_local, yTriple_local]);
                    iE_iTwin_dT_list(:,3) = iE_iTwin_dT_list(:,3) * umPerX;   % --> because unit is 'x/y coordinate' rather than 'data point/index', do not need to use the factor of 5.
                    iiE_each_twin = find_initial_iE_of_twin_in_grain(struCell, ID_current);
                    % repeat for the neighbor
                    [iE_iTwin_dT_list_nb, valid_grain_b] = find_activity_at_boundary_from_struCell(struCell, iE, ID_neighbor, gb, [xTriple_local, yTriple_local]);
                    iE_iTwin_dT_list_nb(:,3) = iE_iTwin_dT_list_nb(:,3) * umPerX;
                    iiE_of_each_twin_nb = find_initial_iE_of_twin_in_grain(struCell, ID_neighbor);
                    
                    exz_iTwin_jMode = calculate_exz(euler, euler_nb,sampleMaterial);
                    exzr = exz_iTwin_jMode./max(exz_iTwin_jMode,[],1);    % because it is normalized to 1, it might be reasonable to use the value (rather than rank) to represent the easiness.
                
                    if (valid_grain_a)&&(valid_grain_b)
                        % (1.2) Find which of the active twin systems in the grain are 'incoming' to this grain boundary
                        % In addition, we need to find out if this twin first appear at this grain boundary.
                        incoming_TS_TF = zeros(1,6);           % IF the twin system is active on the grain of interest side of the boundary
                        iiE_each_iTwin_at_this_boundary = inf * ones(1,6);    % the first iE that this twin shown at this gb
                        an_initiating_boundary_of_twin = zeros(1,6);    % IF the twin was at this boundary at the iE that it just activated, 'an' initiating rather than 'the' initiating
                        intersection_to_triple = nan*zeros(6,1);    % Distance of intersection to triple point
                        twin_at_triple = zeros(6,1);        % In fact no need to label.  Just compare intersection_to_triple to a distance should be enough.  % Does this twin only intersect gb close to triple point
                        
                        for iTwin = 1:6
                            ind = (iE_iTwin_dT_list(:,2)==iTwin)&(iE_iTwin_dT_list(:,3)>dToTriple_th); % This considers TriplePoint, eliminating intersections too close to triple points
                            if any(ind)
                                incoming_TS_TF(iTwin) = 1;
                                iiE_each_iTwin_at_this_boundary(iTwin) = min(iE_iTwin_dT_list(ind,1));
                                if iiE_each_iTwin_at_this_boundary(iTwin) == iiE_each_twin(iTwin)
                                    an_initiating_boundary_of_twin(iTwin) = 1;
                                end
                                
                                dmax = max(iE_iTwin_dT_list(ind,3),[],1);                  % the largest among distance_to_the_closest_triple_point
                                intersection_to_triple(iTwin) = dmax;
                                
                                % If this is still very small, say less than [3 micron] to triple point
                                if dmax < dToTriple_th_to_label
                                    twin_at_triple(iTwin) = 1;
                                    if plotTF
                                        indr = find(dmax == iE_iTwin_dT_list(:,3));
                                        text(iE_iTwin_dT_list(indr,4), iE_iTwin_dT_list(indr,5), 50, [num2str(dmax)], 'color','b','fontsize',36);
                                    end
                                end
                            end
                        end
                        
                        if sum(incoming_TS_TF)==0
                            have_incoming_twin_at_boundary = 0;
                        else
                            have_incoming_twin_at_boundary = 1;
                        end
                        
                        
                        % (1.3) Similarly, find if this neighbor grain (ID_neighbor) is twinned.
                        % If twinned, and twin intersect this_uniqueGB, consider only these intersecting twins
                        % else, consider only basal in this_neighbor
                        % Repeat for the neighbor.  The only difference is that, do not care if it is 'an' initiating boundary.  Just care if it is twinned
                        outgoing_TS_TF = zeros(1,6);           % IF the twin system is active on the grain of interest side of the boundary
                        iiE_each_twin_at_this_boundary_nb = inf * ones(1,6);    % the first iE that this twin shown at this gb
                        intersection_to_triple_nb = nan*zeros(6,1);    % Distance of intersection to triple point
                        twin_at_triple_nb = zeros(6,1);        % Does this twin only intersect gb close to triple point
                        
                        for iTwin = 1:6
                            ind = (iE_iTwin_dT_list_nb(:,2)==iTwin)&(iE_iTwin_dT_list_nb(:,3)>dToTriple_th);	% This considers TriplePoint, eliminating intersections too close to triple points
                            if any(ind)
                                outgoing_TS_TF(iTwin) = 1;     % With d_to_triple considered
                                iiE_each_twin_at_this_boundary_nb(iTwin) = min(iE_iTwin_dT_list_nb(ind,1));
                                
                                dmax = max(iE_iTwin_dT_list_nb(ind,3),[],1);                  % the largest among distance_to_the_closest_triple_point
                                intersection_to_triple_nb(iTwin) = dmax;
                                if dmax < dToTriple_th_to_label
                                    twin_at_triple_nb(iTwin) = 1;
                                    if plotTF==1
                                        indr = find(dmax == iE_iTwin_dT_list_nb(:,3));
                                        text(iE_iTwin_dT_list_nb(indr,4), iE_iTwin_dT_list_nb(indr,5), 50, [num2str(dmax)], 'color','b','fontsize',36);
                                    end
                                end
                            end
                        end
                        
                        
                        if sum(outgoing_TS_TF)==0
                            % Only basal in this_nb
                            outgoing_basal = 1;
                            outgoing_ssn_0 = [1 2 3];
                        else
                            outgoing_basal = 0;
                            tempV = 19:24;
                            outgoing_ssn_0 = tempV(logical(outgoing_TS_TF));
                        end
                        
                        % Calculate mPrime
                        ind_nb = find(gID==ID_neighbor);
                        euler_nb = [gPhi1(ind_nb),gPhi(ind_nb),gPhi2(ind_nb)];
                        [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler, euler_nb, stressTensor, gb_normal_v, sampleMaterial, 'twin');
                        
                        [max_basal_SF, ~] = max(schmidFactorG1(1:3));
                        [max_twin_SF, ~] = max(schmidFactorG1(19:24));
                        [max_basal_SF_nb, ~] = max(schmidFactorG2(1:3));
                        [max_twin_SF_nb, ~] = max(schmidFactorG2(19:24));
                        
                        % Based on whether considering basal or twin in this neighbor, do different things: ...
                        % (3) label basal/twin slip trace
                        ID_target = ID_neighbor;
                        if plotTF==1
                            local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, outgoing_ssn_0, gca);
                        end
                        
                        % find SFs in the neighbor grain
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin,~] = sort(mPrime_each_twin,'descend');
                        [~,rank_mPrime_each_twin] = ismember(mPrime_each_twin, sorted_mPrime_each_twin);    % m'-rank
                        ssn_nb_for_each_twin = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin,~] = sort(resB_each_twin,'ascend');
                        [~,rank_resB_each_twin] = ismember(resB_each_twin, sorted_resB_each_twin);  % resB-rank
                        ssn_nb_for_each_twin_r = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        % The above considers what ss/ts were found experimentally.  
                        % The following considers theoretically, [m', m'-rank, corresponding SF, SF number in neighbor]  
                        % (1) for basals in neighbor  
                        outgoing_ssn_0 = [1 2 3];
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin_wrtB, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin_wrtB,~] = sort(mPrime_each_twin_wrtB,'descend');
                        [~,rank_mPrime_each_twin_wrtB] = ismember(mPrime_each_twin_wrtB, sorted_mPrime_each_twin_wrtB);    % m'-rank
                        ssn_nb_for_each_twin_wrtB = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_wrtB = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin_wrtB, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin_wrtB,~] = sort(resB_each_twin_wrtB,'ascend');
                        [~,rank_resB_each_twin_wrtB] = ismember(resB_each_twin_wrtB, sorted_resB_each_twin_wrtB);  % resB-rank
                        ssn_nb_for_each_twin_r_wrtB = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r_wrtB = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        % (2) for ETWs in neighbor  
                        outgoing_ssn_0 = [19 20 21 22 23 24];
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin_wrtT, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin_wrtT,~] = sort(mPrime_each_twin_wrtT,'descend');
                        [~,rank_mPrime_each_twin_wrtT] = ismember(mPrime_each_twin_wrtT, sorted_mPrime_each_twin_wrtT);    % m'-rank
                        ssn_nb_for_each_twin_wrtT = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_wrtT = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin_wrtT, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin_wrtT,~] = sort(resB_each_twin_wrtT,'ascend');
                        [~,rank_resB_each_twin_wrtT] = ismember(resB_each_twin_wrtT, sorted_resB_each_twin_wrtT);  % resB-rank
                        ssn_nb_for_each_twin_r_wrtT = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r_wrtT = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        
                        
                        % [NOTE] --> Could draw mPrime and color, but label ResB, to maximize information provided.
                        if plotTF==1
                            % find the mPrime to plot. ------> Note that for each boundary, only one mPrime and rank is selected to plot.  But for each twin system, there is an mPrime.
                            if have_incoming_twin_at_boundary
                                tv = double(incoming_TS_TF(:)); % temp vector
                                tv(tv==0)=nan;
                                [mPrime_toPlot, ind] = max(mPrime_each_twin.*tv);
                                mPrime_rank_toPlot = rank_mPrime_each_twin(ind);
                            else
                                tv = double(activeTS(:));
                                tv(tv==0)=nan;
                                [mPrime_toPlot, ind] = max(mPrime_each_twin.*tv);
                                mPrime_rank_toPlot = rank_mPrime_each_twin(ind);
                            end
                            
                            % find the resB to plot
                            if have_incoming_twin_at_boundary
                                [resB_toPlot, ind] = min(resB_each_twin./incoming_TS_TF(:));
                                resB_rank_toPlot = rank_resB_each_twin(ind);
                            else
                                [resB_toPlot, ind] = min(resB_each_twin./activeTS(:));
                                resB_rank_toPlot = rank_resB_each_twin(ind);
                            end
                            
                            
                            inds = find(uniqueBoundary_local==gb);
                            plot(X_local(inds),Y_local(inds),'.','markersize',18,'color',[1 1 0] * (mPrime_toPlot+1)/2);
                            if outgoing_basal==1
                                text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','r','fontsize',16);
                            else
                                text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','b','fontsize',16);
                            end
                        end
                        

                        % Assign activity type, and calculate strain of interest
                        % need to determine 
                        % iE_g: the earliest iE that a twin was active at this gb in the grain of interest, and
                        % iE_n: the earliest iE that a twin was active at this gb in the neighbor of interest 
                        ir = iE_iTwin_dT_list(:,3) > dToTriple_th;
                        iE_g = min(iE_iTwin_dT_list(ir,1));
                        ir = iE_iTwin_dT_list_nb(:,3) > dToTriple_th;
                        iE_n = min(iE_iTwin_dT_list_nb(ir,1));
                        % Artificially make it 'inf' in case it is empty  
                        if isempty(iE_g)
                            iE_g = inf;
                        end
                        if isempty(iE_n)
                            iE_n = inf;
                        end

                        
                        if (iE_g>iE)%||(isempty(iE_g))
                            if (iE_n>iE)%||(isempty(iE_n))
                                % not involved
                                bg_not_involved = [bg_not_involved; gb, ID_current];
                            elseif (iE_n==iE)
                                % slip induced twin, slip side
                                bg_slip_twin_a = [bg_slip_twin_a; gb, ID_current];
                            elseif (iE_n<iE)
                                % slip induced twin growth, slip side
                                bg_slip_growth_a = [bg_slip_growth_a; gb, ID_current];
                            end
                        elseif (iE_g==iE)
                            % just twinned this grain
                            if (iE_n>iE)%||(isempty(iE_n))
                                % slip induced twin, twin side
                                bg_slip_twin_b = [bg_slip_twin_b; gb, ID_current];
                            elseif (iE_n==iE)
                                % co-found twin
                                bg_co_found = [bg_co_found; gb, ID_current];
                            elseif (iE_n<iE)
                                % twin induced twin, new twin side
                                bg_twin_twin_b = [bg_twin_twin_b; gb, ID_current];
                            end
                        elseif (iE_g<iE)
                            % twinned this grain
                            if (iE_n>iE)%||(isempty(iE_n))
                                % slip induced twin growth, twin side
                                bg_slip_growth_b = [bg_slip_growth_b; gb, ID_current];
                            elseif (iE_n==iE)
                                % twin induced twin, old twin side
                                bg_twin_twin_a = [bg_twin_twin_a; gb, ID_current];
                            elseif (iE_n<iE)
                                % co-growth
                                bg_co_growth = [bg_co_growth; gb, ID_current];
                            end
                        end
                        
                        
                        % Extend table for summary.  -- can assign in local table using name?
                        for iTwin=1:6
                            lookupMa = [lookupMa; ID_current, ID_neighbor, iTwin, tSF(iTwin), reshape(schmidFactorG2,1,[]), mPrimeMatrix(iTwin+18,:)];
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
                            T_local.gb_length(ir) = gb_length;
                            T_local.gb_dir(ir) = gb_dir;
                            
                            T_local.incoming(ir) = incoming_TS_TF(iTwin);
                            T_local.iiE_each_twin(ir) = iiE_each_twin(iTwin);
                            T_local.iiE_each_twin_at_this_boundary(ir) = iiE_each_iTwin_at_this_boundary(iTwin);
                            T_local.intersection_to_triple(ir) = intersection_to_triple(iTwin);
                            T_local.iiE_twins_at_this_boundary_nb(ir) = min(iiE_each_twin_at_this_boundary_nb);
                            
                            T_local.mPrime(ir) = mPrime_each_twin(iTwin);
                            T_local.rank_mPrime(ir) = rank_mPrime_each_twin(iTwin);
                            T_local.ssn_nb(ir) = ssn_nb_for_each_twin(iTwin);
                            T_local.SF_nb(ir) = SF_nb_for_each_twin(iTwin);
                            
                            T_local.resB(ir) = resB_each_twin(iTwin);
                            T_local.rank_resB(ir) = rank_resB_each_twin(iTwin);
                            T_local.ssn_nb_r(ir) = ssn_nb_for_each_twin_r(iTwin);
                            T_local.SF_nb_r(ir) = SF_nb_for_each_twin_r(iTwin);
                            
                            T_local.mPrime_wrtB(ir) = mPrime_each_twin_wrtB(iTwin);
                            T_local.rank_mPrime_wrtB(ir) = rank_mPrime_each_twin_wrtB(iTwin);
                            T_local.ssn_nb_wrtB(ir) = ssn_nb_for_each_twin_wrtB(iTwin);
                            T_local.SF_nb_wrtB(ir) = SF_nb_for_each_twin_wrtB(iTwin);
                            
                            T_local.resB_wrtB(ir) = resB_each_twin_wrtB(iTwin);
                            T_local.rank_resB_wrtB(ir) = rank_resB_each_twin_wrtB(iTwin);
                            T_local.ssn_nb_r_wrtB(ir) = ssn_nb_for_each_twin_r_wrtB(iTwin);
                            T_local.SF_nb_r_wrtB(ir) = SF_nb_for_each_twin_r_wrtB(iTwin);
                            
                            T_local.mPrime_wrtT(ir) = mPrime_each_twin_wrtT(iTwin);
                            T_local.rank_mPrime_wrtT(ir) = rank_mPrime_each_twin_wrtT(iTwin);
                            T_local.ssn_nb_wrtT(ir) = ssn_nb_for_each_twin_wrtT(iTwin);
                            T_local.SF_nb_wrtT(ir) = SF_nb_for_each_twin_wrtT(iTwin);
                            
                            T_local.resB_wrtT(ir) = resB_each_twin_wrtT(iTwin);
                            T_local.rank_resB_wrtT(ir) = rank_resB_each_twin_wrtT(iTwin);
                            T_local.ssn_nb_r_wrtT(ir) = ssn_nb_for_each_twin_r_wrtT(iTwin);
                            T_local.SF_nb_r_wrtT(ir) = SF_nb_for_each_twin_r_wrtT(iTwin);
                            
                            T_local.initiating(ir) = false;
                            
                            T_local.eMean_1(ir) = eMean_1;
                            T_local.eMean_2(ir) = eMean_2;
                            T_local.eMean_1_nb(ir) = eMean_1_nb;
                            T_local.eMean_2_nb(ir) = eMean_2_nb;
                            T_local.max_basal_SF(ir) = max_basal_SF;
                            T_local.max_twin_SF(ir) = max_twin_SF;
                            T_local.max_basal_SF_nb(ir) = max_basal_SF_nb;
                            T_local.max_twin_SF_nb(ir) = max_twin_SF_nb;
                            
                            T_local.exz_ba(ir) = exz_iTwin_jMode(iTwin, 1);
                            T_local.exz_pr(ir) = exz_iTwin_jMode(iTwin, 2);
                            T_local.exz_py(ir) = exz_iTwin_jMode(iTwin, 3);
                            T_local.exz_pyII(ir) = exz_iTwin_jMode(iTwin, 4);
                            T_local.exz_etw(ir) = exz_iTwin_jMode(iTwin, 5);
                            
                            T_local.exzr_ba(ir) = exzr(iTwin, 1);
                            T_local.exzr_pr(ir) = exzr(iTwin, 2);
                            T_local.exzr_py(ir) = exzr(iTwin, 3);
                            T_local.exzr_pyII(ir) = exzr(iTwin, 4);
                            T_local.exzr_etw(ir) = exzr(iTwin, 5);
                            
                            ind = find(struCell{iE}(iS).tGb{iTwin} == gb);
                            if isempty(ind)
                                tGbVol = 0;
                            else
                                tGbVol = struCell{iE}(iS).tGbVol{iTwin}(ind);
                            end
                            
                            % add to new colume
                            T_local.tGbVol(ir) = tGbVol;
                            T_local.tGbVolPct(ir) = tGbVol/(pi/4*(gDia/umPerDp)^2);
                            T_local.tGbStrength(ir) = T_local.tGbVolPct(ir)/gb_length;
                        end

                    else
                        disp('not valid_grain_a and valid_grain_b')
                    end
                    % end of (valid_grain_a)&&(valid_grain_b)
                    
                end
                % end of ~isempty(iS_neighbor)
                
            end
            % end of for iNb = 1:nNeighbors
            
            % 'incoming_TS' just means intersecting.  Each twin can intersect several gbs.
            % We want to find the one with the (highest m') to represent the initiating point.
            % But there are potential pairs of slip/twin we can choose from, e.g., with a threshold for Schmid factor, and grain size, etc  
            for iTwin = 1:6
                rows = find((T_local.TS==iTwin)&(T_local.iiE_each_twin==T_local.iiE_each_twin_at_this_boundary)&(T_local.gDia>80)&(T_local.gDia_neighbor>80));
                if ~isempty(rows)
                    [~,ind] = max(T_local.mPrime(rows,:));
                    T_local.initiating(rows(ind)) = 1;
                    % disp(['Twin ',num2str(iTwin),' first appear, initiatiated from neighbor ID:',num2str(T_local.ID_neighbor(rows(ind)))]);
                end
            end
            ind = (T_local.initiating==1);
%            label_map_with_ID(X_local,Y_local,ID_local,gcf,T_local.ID_neighbor(ind),'m');
            T = [T;T_local];
            
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
            
        else
%             % can use continue to just plot twinned ones
%             iS = iS + 1;
%             continue;
            
            % If grain of interest is not twinned, this means (iE_g > iE)   
            
            T_local = T_template;
            % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
            for iNb = 1:nNeighbors
                ID_neighbor = ID_neighbors(iNb);
                iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
                
                ind = find(gID==ID_neighbor);
                euler_nb = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
                
                if ~isempty(iS_neighbor)
                    gDia_neighbor = sqrt(4*(struCell{iE}(iS_neighbor).gVol* (360/4096*5)^2)/pi);
                    % activeTS_nb = sum(struCell{iE}(iS_neighbor).cTrueTwin,1)>0;
                    activeTS_nb = cellfun(@(x) ~isempty(x), struCell{iE}(iS_neighbor).tGb);
                    
                    
                    % (1.1) Calculate this_uniqueGB number.
                    if ID_current > ID_neighbor
                        gb = ID_current * 10000 + ID_neighbor;
                    else
                        gb = ID_neighbor * 10000 + ID_current;
                    end
                    % Find gb_length and gb_dir
                    ind = uniqueBoundary_local(:)==gb;
                    gb_length = sum(ind);
                    xs = X_local(ind);
                    ys = Y_local(ind);
                    mdl = fitlm(xs,ys);
                    gb_dir = atand(mdl.Coefficients.Estimate(2));
                    gb_normal = 90 + gb_dir;
                    gb_normal_v = [cosd(gb_normal),sind(gb_normal),0]; 
                    gb_normal_v = gb_normal_v/norm(gb_normal_v);
                    
                    % strain calculation in area of interest.
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                    ind = (distMap_local>0)&(distMap_local<=50);
                    eMean_1 =  nanmean(eMap_local(ind));   % mean
                    ind = (distMap_local>50)&(distMap_local<=100);
                    eMean_2 =  nanmean(eMap_local(ind));   % mean
                    
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_neighbor]);
                    ind = (distMap_local>0)&(distMap_local<=50);
                    eMean_1_nb =  nanmean(eMap_local(ind));   % mean
                    ind = (distMap_local>50)&(distMap_local<=100);
                    eMean_2_nb =  nanmean(eMap_local(ind));   % mean
                    
                    ind_of_indTriple = sum(ismember(tripleIDs,ID_current)+ismember(tripleIDs,ID_neighbor), 2)==2;
                    xTriple_local = X(indTriple(ind_of_indTriple));
                    yTriple_local = Y(indTriple(ind_of_indTriple));
                    
                    % Find [iE, iTwin, dToTriple] on the grain_of_interest side, then on the neighboring side
                    [iE_iTwin_dT_list, valid_grain_a] = find_activity_at_boundary_from_struCell(struCell, iE, ID_current, gb, [xTriple_local, yTriple_local]);
                    iE_iTwin_dT_list(:,3) = iE_iTwin_dT_list(:,3) * umPerX;   % --> because unit is 'x/y coordinate' rather than 'data point/index', do not need to use the factor of 5.
                    iiE_each_twin = find_initial_iE_of_twin_in_grain(struCell, ID_current);
                    % repeat for the neighbor
                    [iE_iTwin_dT_list_nb, valid_grain_b] = find_activity_at_boundary_from_struCell(struCell, iE, ID_neighbor, gb, [xTriple_local, yTriple_local]);
                    iE_iTwin_dT_list_nb(:,3) = iE_iTwin_dT_list_nb(:,3) * umPerX;
                    iiE_of_each_twin_nb = find_initial_iE_of_twin_in_grain(struCell, ID_neighbor);
                    
                    exz_iTwin_jMode = calculate_exz(euler, euler_nb, sampleMaterial);
                    exzr = exz_iTwin_jMode./max(exz_iTwin_jMode,[],1);    % because it is normalized to 1, it might be reasonable to use the value (rather than rank) to represent the easiness.
                
                    if (valid_grain_a)&&(valid_grain_b)
                        % (1.2) Find which of the active twin systems in the grain are 'incoming' to this grain boundary
                        % In addition, we need to find out if this twin first appear at this grain boundary.
                        incoming_TS_TF = zeros(1,6);           % IF the twin system is active on the grain of interest side of the boundary
                        iiE_each_iTwin_at_this_boundary = inf * ones(1,6);    % the first iE that this twin shown at this gb
                        an_initiating_boundary_of_twin = zeros(1,6);    % IF the twin was at this boundary at the iE that it just activated, 'an' initiating rather than 'the' initiating
                        intersection_to_triple = nan*zeros(6,1);    % Distance of intersection to triple point
                        twin_at_triple = zeros(6,1);        % In fact no need to label.  Just compare intersection_to_triple to a distance should be enough.  % Does this twin only intersect gb close to triple point
                        
                        % no twin is active in parent grain
                        for iTwin = 1:6
                            ind = (iE_iTwin_dT_list(:,2)==iTwin)&(iE_iTwin_dT_list(:,3)>dToTriple_th); % This considers TriplePoint, eliminating intersections too close to triple points
                            if any(ind)
                                error('should not have any twin system in grain of interest here'); 
                                incoming_TS_TF(iTwin) = 1;
                                iiE_each_iTwin_at_this_boundary(iTwin) = min(iE_iTwin_dT_list(ind,1));
                                if iiE_each_iTwin_at_this_boundary(iTwin) == iiE_each_twin(iTwin)
                                    an_initiating_boundary_of_twin(iTwin) = 1;
                                end
                                
                                dmax = max(iE_iTwin_dT_list(ind,3),[],1);                  % the largest among distance_to_the_closest_triple_point
                                intersection_to_triple(iTwin) = dmax;
                                
                                % If this is still very small, say less than [3 micron] to triple point
                                if dmax < dToTriple_th_to_label
                                    twin_at_triple(iTwin) = 1;
                                    if plotTF
                                        indr = find(dmax == iE_iTwin_dT_list(:,3));
                                        text(iE_iTwin_dT_list(indr,4), iE_iTwin_dT_list(indr,5), 50, [num2str(dmax)], 'color','b','fontsize',36);
                                    end
                                end
                            end
                        end
                        
                        if sum(incoming_TS_TF)==0
                            have_incoming_twin_at_boundary = 0;
                        else
                            have_incoming_twin_at_boundary = 1;
                        end
                        
                        
                        % (1.3) Similarly, find if this neighbor grain (ID_neighbor) is twinned.
                        % If twinned, and twin intersect this_uniqueGB, consider only these intersecting twins
                        % else, consider only basal in this_neighbor
                        % Repeat for the neighbor.  The only difference is that, do not care if it is 'an' initiating boundary.  Just care if it is twinned
                        outgoing_TS_TF = zeros(1,6);           % IF the twin system is active on the grain of interest side of the boundary
                        iiE_each_twin_at_this_boundary_nb = inf * ones(1,6);    % the first iE that this twin shown at this gb
                        intersection_to_triple_nb = nan*zeros(6,1);    % Distance of intersection to triple point
                        twin_at_triple_nb = zeros(6,1);        % Does this twin only intersect gb close to triple point
                        
                        for iTwin = 1:6
                            ind = (iE_iTwin_dT_list_nb(:,2)==iTwin)&(iE_iTwin_dT_list_nb(:,3)>dToTriple_th);	% This considers TriplePoint, eliminating intersections too close to triple points
                            if any(ind)
                                outgoing_TS_TF(iTwin) = 1;     % With d_to_triple considered
                                iiE_each_twin_at_this_boundary_nb(iTwin) = min(iE_iTwin_dT_list_nb(ind,1));
                                
                                dmax = max(iE_iTwin_dT_list_nb(ind,3),[],1);                  % the largest among distance_to_the_closest_triple_point
                                intersection_to_triple_nb(iTwin) = dmax;
                                if dmax < dToTriple_th_to_label
                                    twin_at_triple_nb(iTwin) = 1;
%                                     if plotTF==1
%                                         indr = find(dmax == iE_iTwin_dT_list_nb(:,3));
%                                         text(iE_iTwin_dT_list_nb(indr,4), iE_iTwin_dT_list_nb(indr,5), 50, [num2str(dmax)], 'color','b','fontsize',36);
%                                     end
                                end
                            end
                        end
                        
                        
                        if sum(outgoing_TS_TF)==0
                            % Only basal in this_nb
                            outgoing_basal = 1;
                            outgoing_ssn_0 = [1 2 3];
                        else
                            outgoing_basal = 0;
                            tempV = 19:24;
                            outgoing_ssn_0 = tempV(logical(outgoing_TS_TF));
                        end
                        
                        % Calculate mPrime
                        ind_nb = find(gID==ID_neighbor);
                        euler_nb = [gPhi1(ind_nb),gPhi(ind_nb),gPhi2(ind_nb)];
                        [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler, euler_nb, stressTensor, gb_normal_v, sampleMaterial, 'twin');
                        
                        [max_basal_SF, ~] = max(schmidFactorG1(1:3));
                        [max_twin_SF, ~] = max(schmidFactorG1(19:24));
                        [max_basal_SF_nb, ~] = max(schmidFactorG2(1:3));
                        [max_twin_SF_nb, ~] = max(schmidFactorG2(19:24));
                        
                        % Based on whether considering basal or twin in this neighbor, do different things: ...
                        % (3) label basal/twin slip trace
                        ID_target = ID_neighbor;
%                         if plotTF==1
%                             local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, outgoing_ssn_0, gca);
%                         end
                        
                        % find SFs in the neighbor grain
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.  
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin,~] = sort(mPrime_each_twin,'descend');
                        [~,rank_mPrime_each_twin] = ismember(mPrime_each_twin, sorted_mPrime_each_twin);    % m'-rank
                        ssn_nb_for_each_twin = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin,~] = sort(resB_each_twin,'ascend');
                        [~,rank_resB_each_twin] = ismember(resB_each_twin, sorted_resB_each_twin);  % resB-rank
                        ssn_nb_for_each_twin_r = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        % The above considers what ss/ts were found experimentally.  
                        % The following considers theoretically, [m', m'-rank, corresponding SF, SF number in neighbor]  
                        % (1) for basals in neighbor  
                        outgoing_ssn_0 = [1 2 3];
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin_wrtB, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin_wrtB,~] = sort(mPrime_each_twin_wrtB,'descend');
                        [~,rank_mPrime_each_twin_wrtB] = ismember(mPrime_each_twin_wrtB, sorted_mPrime_each_twin_wrtB);    % m'-rank
                        ssn_nb_for_each_twin_wrtB = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_wrtB = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin_wrtB, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin_wrtB,~] = sort(resB_each_twin_wrtB,'ascend');
                        [~,rank_resB_each_twin_wrtB] = ismember(resB_each_twin_wrtB, sorted_resB_each_twin_wrtB);  % resB-rank
                        ssn_nb_for_each_twin_r_wrtB = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r_wrtB = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        % (2) for ETWs in neighbor  
                        outgoing_ssn_0 = [19 20 21 22 23 24];
                        SFs_nb_0 = schmidFactorG2(outgoing_ssn_0);
                        % In addition, not all ssn_nb are allowed to analyze its m' value
                        % Maybe it's better to find potential systems, e.g., all those with sf > [0.2], or normalized > 0.8
                        SF_nb_enough_TF = (SFs_nb_0>0.2) | ((SFs_nb_0>0)&(SFs_nb_0./max(SFs_nb_0)>0.8));   % logical vector, indicating wheter the out-going ss is considered as able to be activated (e.g., SF high enough)
                        % [detail] sometimes all basal can have SF=0.
                        if sum(SF_nb_enough_TF)==0
                            SF_nb_enough_TF(:) = true;
                        end
                        outgoing_ssn = outgoing_ssn_0(SF_nb_enough_TF);
                        SFs_nb = schmidFactorG2(outgoing_ssn);
                        
                        mPrime_local = mPrimeMatrix(19:24, outgoing_ssn);
                        [mPrime_each_twin_wrtT, ind] = max(mPrime_local,[],2);     % value of interests: m'
                        [sorted_mPrime_each_twin_wrtT,~] = sort(mPrime_each_twin_wrtT,'descend');
                        [~,rank_mPrime_each_twin_wrtT] = ismember(mPrime_each_twin_wrtT, sorted_mPrime_each_twin_wrtT);    % m'-rank
                        ssn_nb_for_each_twin_wrtT = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_wrtT = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest

                        resB_local = resBurgersMatrix(19:24, outgoing_ssn);
                        [resB_each_twin_wrtT, ind] = min(resB_local,[],2);       % value of interest: resB
                        [sorted_resB_each_twin_wrtT,~] = sort(resB_each_twin_wrtT,'ascend');
                        [~,rank_resB_each_twin_wrtT] = ismember(resB_each_twin_wrtT, sorted_resB_each_twin_wrtT);  % resB-rank
                        ssn_nb_for_each_twin_r_wrtT = outgoing_ssn(ind); % ssn in neighbor, for each ss/ts in the grain of interest
                        SF_nb_for_each_twin_r_wrtT = SFs_nb(ind);        % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                        
                        
                        
                        % [NOTE] --> Could draw mPrime and color, but label ResB, to maximize information provided.
%                         if plotTF==1
%                             % find the mPrime to plot. ------> Note that for each boundary, only one mPrime and rank is selected to plot.  But for each twin system, there is an mPrime.
%                             if have_incoming_twin_at_boundary
%                                 [mPrime_toPlot, ind] = max(mPrime_each_twin.*incoming_TS_TF(:));
%                                 mPrime_rank_toPlot = rank_mPrime_each_twin(ind);
%                             else
%                                 [mPrime_toPlot, ind] = max(mPrime_each_twin.*activeTS(:));
%                                 mPrime_rank_toPlot = rank_mPrime_each_twin(ind);
%                             end
%                             
%                             % find the resB to plot
%                             if have_incoming_twin_at_boundary
%                                 [resB_toPlot, ind] = min(resB_each_twin./incoming_TS_TF(:));
%                                 resB_rank_toPlot = rank_resB_each_twin(ind);
%                             else
%                                 [resB_toPlot, ind] = min(resB_each_twin./activeTS(:));
%                                 resB_rank_toPlot = rank_resB_each_twin(ind);
%                             end
%                             
%                             
%                             inds = find(uniqueBoundary_local==gb);
%                             plot(X_local(inds),Y_local(inds),'.','markersize',18,'color',[1 1 0] * mPrime_toPlot);
%                             if outgoing_basal==1
%                                 text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_rank_toPlot + mPrime_toPlot,4)], 'color','r','fontsize',16);
%                             else
%                                 text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_rank_toPlot + mPrime_toPlot,4)], 'color','b','fontsize',16);
%                             end
%                         end
                        

                        % iE_g > iE
                        ir = iE_iTwin_dT_list_nb(:,3) > dToTriple_th;
                        iE_n = min(iE_iTwin_dT_list_nb(ir,1));
                        if isempty(iE_n)
                            iE_n = inf;
                        end

                        if (iE_n>iE)%||(isempty(iE_n))
                            % not involved
                            bg_not_involved = [bg_not_involved; gb, ID_current];
                        elseif (iE_n==iE)
                            % slip induced twin, slip side
                            bg_slip_twin_a = [bg_slip_twin_a; gb, ID_current];
                        elseif (iE_n<iE)
                            % slip induced twin growth, slip side
                            bg_slip_growth_a = [bg_slip_growth_a; gb, ID_current];
                        end

                        
                        
                        % Extend table for summary.  -- can assign in local table using name?
                        for iTwin=1:6
                            lookupMb = [lookupMb; ID_current, ID_neighbor, iTwin, tSF(iTwin), reshape(schmidFactorG2,1,[]), mPrimeMatrix(iTwin+18,:)];
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
                            T_local.gb_length(ir) = gb_length;
                            T_local.gb_dir(ir) = gb_dir;
                                                        
                            T_local.incoming(ir) = incoming_TS_TF(iTwin);
                            T_local.iiE_each_twin(ir) = iiE_each_twin(iTwin);
                            T_local.iiE_each_twin_at_this_boundary(ir) = iiE_each_iTwin_at_this_boundary(iTwin);
                            T_local.intersection_to_triple(ir) = intersection_to_triple(iTwin);
                            T_local.iiE_twins_at_this_boundary_nb(ir) = min(iiE_each_twin_at_this_boundary_nb);
                            
                            T_local.mPrime(ir) = mPrime_each_twin(iTwin);
                            T_local.rank_mPrime(ir) = rank_mPrime_each_twin(iTwin);
                            T_local.ssn_nb(ir) = ssn_nb_for_each_twin(iTwin);
                            T_local.SF_nb(ir) = SF_nb_for_each_twin(iTwin);
                            
                            T_local.resB(ir) = resB_each_twin(iTwin);
                            T_local.rank_resB(ir) = rank_resB_each_twin(iTwin);
                            T_local.ssn_nb_r(ir) = ssn_nb_for_each_twin_r(iTwin);
                            T_local.SF_nb_r(ir) = SF_nb_for_each_twin_r(iTwin);
                            
                            T_local.mPrime_wrtB(ir) = mPrime_each_twin_wrtB(iTwin);
                            T_local.rank_mPrime_wrtB(ir) = rank_mPrime_each_twin_wrtB(iTwin);
                            T_local.ssn_nb_wrtB(ir) = ssn_nb_for_each_twin_wrtB(iTwin);
                            T_local.SF_nb_wrtB(ir) = SF_nb_for_each_twin_wrtB(iTwin);
                            
                            T_local.resB_wrtB(ir) = resB_each_twin_wrtB(iTwin);
                            T_local.rank_resB_wrtB(ir) = rank_resB_each_twin_wrtB(iTwin);
                            T_local.ssn_nb_r_wrtB(ir) = ssn_nb_for_each_twin_r_wrtB(iTwin);
                            T_local.SF_nb_r_wrtB(ir) = SF_nb_for_each_twin_r_wrtB(iTwin);
                            
                            T_local.mPrime_wrtT(ir) = mPrime_each_twin_wrtT(iTwin);
                            T_local.rank_mPrime_wrtT(ir) = rank_mPrime_each_twin_wrtT(iTwin);
                            T_local.ssn_nb_wrtT(ir) = ssn_nb_for_each_twin_wrtT(iTwin);
                            T_local.SF_nb_wrtT(ir) = SF_nb_for_each_twin_wrtT(iTwin);
                            
                            T_local.resB_wrtT(ir) = resB_each_twin_wrtT(iTwin);
                            T_local.rank_resB_wrtT(ir) = rank_resB_each_twin_wrtT(iTwin);
                            T_local.ssn_nb_r_wrtT(ir) = ssn_nb_for_each_twin_r_wrtT(iTwin);
                            T_local.SF_nb_r_wrtT(ir) = SF_nb_for_each_twin_r_wrtT(iTwin);
                            
                            T_local.initiating(ir) = false;
                            
                            T_local.eMean_1(ir) = eMean_1;
                            T_local.eMean_2(ir) = eMean_2;
                            T_local.eMean_1_nb(ir) = eMean_1_nb;
                            T_local.eMean_2_nb(ir) = eMean_2_nb;
                            T_local.max_basal_SF(ir) = max_basal_SF;
                            T_local.max_twin_SF(ir) = max_twin_SF;
                            T_local.max_basal_SF_nb(ir) = max_basal_SF_nb;
                            T_local.max_twin_SF_nb(ir) = max_twin_SF_nb;
                            
                            T_local.exz_ba(ir) = exz_iTwin_jMode(iTwin, 1);
                            T_local.exz_pr(ir) = exz_iTwin_jMode(iTwin, 2);
                            T_local.exz_py(ir) = exz_iTwin_jMode(iTwin, 3);
                            T_local.exz_pyII(ir) = exz_iTwin_jMode(iTwin, 4);
                            T_local.exz_etw(ir) = exz_iTwin_jMode(iTwin, 5);
                            
                            T_local.exzr_ba(ir) = exzr(iTwin, 1);
                            T_local.exzr_pr(ir) = exzr(iTwin, 2);
                            T_local.exzr_py(ir) = exzr(iTwin, 3);
                            T_local.exzr_pyII(ir) = exzr(iTwin, 4);
                            T_local.exzr_etw(ir) = exzr(iTwin, 5);
                            
                            ind = find(struCell{iE}(iS).tGb{iTwin} == gb);
                            if isempty(ind)
                                tGbVol = 0;
                            else
                                tGbVol = struCell{iE}(iS).tGbVol{iTwin}(ind);
                            end
                            
                            % add to new colume
                            T_local.tGbVol(ir) = tGbVol;
                            T_local.tGbVolPct(ir) = tGbVol/(pi/4*(gDia/umPerDp)^2);
                            T_local.tGbStrength(ir) = T_local.tGbVolPct(ir)/gb_length;
                        end

                    else
                        disp('not valid_grain_a and valid_grain_b')
                    end
                    % end of (valid_grain_a)&&(valid_grain_b)
                    
                end
                % end of ~isempty(iS_neighbor)
                
            end
            % end of for iNb = 1:nNeighbors
            
            % 'incoming_TS' just means intersecting.  Each twin can intersect several gbs.
            % We want to find the one with the (highest m') to represent the initiating point.
            % But there are potential pairs of slip/twin we can choose from, e.g., with a threshold for Schmid factor, and grain size, etc  
            for iTwin = 1:6
                rows = find((T_local.TS==iTwin)&(T_local.iiE_each_twin==T_local.iiE_each_twin_at_this_boundary)&(T_local.gDia>80)&(T_local.gDia_neighbor>80));
                if ~isempty(rows)
                    error('should not find initiating twin here.');
                    [~,ind] = max(T_local.mPrime(rows,:));
                    T_local.initiating(rows(ind)) = 1;
                    % disp(['Twin ',num2str(iTwin),' first appear, initiatiated from neighbor ID:',num2str(T_local.ID_neighbor(rows(ind)))]);
                end
            end
            ind = (T_local.initiating==1);
%            label_map_with_ID(X_local,Y_local,ID_local,gcf,T_local.ID_neighbor(ind),'m');
            T2 = [T2; T_local];

            
        end
        
        % disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
        iS = iS + 1;
    end
    close(hW);
    warning on;
    
    %% Save
        timeStr = datestr(now,'yyyymmdd_HHMM');
        save(['temp_results/',timeStr,'_twin_gb_summary_',num2str(iE),'.mat'], 'struCell','T','T2', 'lookupMa','lookupMb', ...
            'bg_not_involved','bg_slip_twin_a','bg_slip_twin_b','bg_co_found','bg_twin_twin_a','bg_twin_twin_b','bg_slip_growth_a','bg_slip_growth_b','bg_co_growth',...
            '-v7.3');
        copyfile(['temp_results/',timeStr,'_twin_gb_summary_',num2str(iE),'.mat'], ['temp_results/twin_gb_summary_',num2str(iE),'.mat']);
    
end

%% To study The effect of high strain on twin growth/activation:*
% For each grain and grain boundary, we calculate the mean of the effective 
% strain of all data points at the same distance (d) to the grain boundary,
% for all distances d=1,2,3,...,250 
% And then for each category, we take all the grain and grain boundaries,
% and take the mean of all the mean effective strains @ each distances to
% grain boundary
% So, we get a curve that represents the strain distribution w.r.t to
% distance to gb, for each category
%% Method-1. Loop each grain, each neighbor (unique gb), calculate distance map. Then summarize e distribution.

% calculate edmat: strain distribution matrix
% [iE, ID_current, gb, nanmean of eEff @ d=1:250 data point distance to that gb.  
x_dist = 1:250;
for iE = 2:5
    
    eMap = calculate_effective_strain(strainFile{iE-1}.exx, strainFile{iE-1}.exy, strainFile{iE-1}.eyy);
    edmat = [];
    
    iS = 1;
    warning('off','MATLAB:table:RowsAddedExistingVars');
    continueTF = true;
    dToTriple_th = 5;       % eliminate intersection whose distance to triple point is smaller than this value
    dToTriple_th_to_label = dToTriple_th;  % label if distance of intersection to triple point is smaller than this value
    
    hW = waitbar(0, ['iE=',num2str(iE),' analyze each grain']);
    hN = length(struCell{iE});
    while (continueTF)&&(iS<=length(struCell{iE}))
        waitbar(iS/hN, hW);
        
        close all;
        ID_current = struCell{iE}(iS).gID
        ind = find(gID==ID_current);
        
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);
        
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);

        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);  % This is for effective strain
        
        % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
        for iNb = 1:nNeighbors
            ID_neighbor = ID_neighbors(iNb);
            iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
            if ~isempty(iS_neighbor)
                
                % (1.1) Calculate this_uniqueGB number.
                if ID_current > ID_neighbor
                    gb = ID_current * 10000 + ID_neighbor;
                else
                    gb = ID_neighbor * 10000 + ID_current;
                end
                
                % strain calculation in area of interest.
                distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];
                
            end
            % end of ~isempty(iS_neighbor)
            
        end

        % disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
        iS = iS + 1;
    end
    close(hW);
    warning on;
    
    switch iE
        case 2
            edmat_2 = edmat;
            save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmat_2','x_dist', '-append','-v7.3');
        case 3
            edmat_3 = edmat;
            save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmat_3','x_dist', '-append','-v7.3');
        case 4
            edmat_4 = edmat;
            save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmat_4','x_dist', '-append','-v7.3');
        case 5
            edmat_5 = edmat;
            save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmat_5','x_dist', '-append','-v7.3');
    end
    
end
%% Method-2.  This methods seems not as good as Method-1, so maybe just for debug, no need to run this part.
% First calculate a distMap. Each data point is affected only by the nearest unique gb.
% Then, Loop each grain, each neighbor (unique gb), crop the distance map. Then summarize e distribution.

% calculate edmat: strain distribution matrix
% [iE, ID_current, gb, nanmean of eEff @ d=1:250 data point distance to that gb.  

% [~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
% uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
% 
% [distMap, gbLabel] = city_block(uniqueBoundary);
% 
% for iE = 2:5
%     
%     eMap = calculate_effective_strain(strainFile{iE-1}.exx, strainFile{iE-1}.exy, strainFile{iE-1}.eyy);
%     edmat = [];
%     
%     iS = 1;
%     warning('off','MATLAB:table:RowsAddedExistingVars');
%     continueTF = true;
%     dToTriple_th = 5;       % eliminate intersection whose distance to triple point is smaller than this value
%     dToTriple_th_to_label = dToTriple_th;  % label if distance of intersection to triple point is smaller than this value
%     
%     hW = waitbar(0, ['iE=',num2str(iE),' analyze each grain']);
%     hN = length(struCell{iE});
%     while (continueTF)&&(iS<=length(struCell{iE}))
%         waitbar(iS/hN, hW);
%         
%         close all;
%         ID_current = struCell{iE}(iS).gID
%         ind = find(gID==ID_current);
%         
%         nNeighbors = gNNeighbors(ind);
%         ID_neighbors = gNeighbors(ind, 1:nNeighbors);
%         
%         ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);
%         
%         % Make it one data point wider on each side
%         indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
%         indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
%         indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
%         indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
%         
%         ID_local = ID(indR_min:indR_max, indC_min:indC_max);
% 
%         eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);  % This is for effective strain
%         gbLabel_local = gbLabel(indR_min:indR_max, indC_min:indC_max);
%         
%         % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
%         for iNb = 1:nNeighbors
%             ID_neighbor = ID_neighbors(iNb);
%             iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
%             if ~isempty(iS_neighbor)
%                 
%                 % (1.1) Calculate this_uniqueGB number.
%                 if ID_current > ID_neighbor
%                     gb = ID_current * 10000 + ID_neighbor;
%                 else
%                     gb = ID_neighbor * 10000 + ID_current;
%                 end
%                 
%                 % strain calculation in area of interest.
%                 distMap_local = distMap(indR_min:indR_max, indC_min:indC_max);
%                 mask = (gbLabel_local==gb)&(ID_local==ID_current);
%                 distMap_local(~mask) = nan;
%                 edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];
%                 
%             end
%             % end of ~isempty(iS_neighbor)
%             
%         end
% 
%         % disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
%         iS = iS + 1;
%     end
%     close(hW);
%     warning on;
%     
%     switch iE
%         case 2
%             edmatII_2 = edmat;
%             save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmatII_2','x_dist', '-append','-v7.3');
%         case 3
%             edmatII_3 = edmat;
%             save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmatII_3','x_dist', '-append','-v7.3');
%         case 4
%             edmatII_4 = edmat;
%             save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmatII_4','x_dist', '-append','-v7.3');
%         case 5
%             edmatII_5 = edmat;
%             save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'edmatII_5','x_dist', '-append','-v7.3');
%     end
%     
% end


%% Added something to the code, but run here to increase speed, no need to calculate everything again.  If run from beginning, no need to run again.
if 0
    for iE = 2:5
        load(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'T','T2');
        for ir = 1:size(T,1)
            ID_current = T.ID(ir)
            ID_neighbor = T.ID_neighbor(ir);
            gbNum = max(ID_current,ID_neighbor)*10000 + min(ID_current,ID_neighbor);
            iTwin = T.TS(ir);
            
            iS = find(arrayfun(@(x) x.gID == ID_current,struCell{iE}));
            
            ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
            if isempty(ind)
                tGbVol = 0;
            else
                tGbVol = struCell{iE}(iS).tGbVol{iTwin}(ind);
            end
            
            % add to new colume
            T.tGbVol(ir) = tGbVol;
            T.tGbVolPct(ir) = tGbVol/(pi/4*(T.gDia(ir)/umPerDp)^2);
            T.tGbStrength(ir) = T.tGbVolPct(ir)/T.gb_length(ir);
        end
        for ir = 1:size(T2,1)
            ID_current = T2.ID(ir)
            ID_neighbor = T2.ID_neighbor(ir);
            gbNum = max(ID_current,ID_neighbor)*10000 + min(ID_current,ID_neighbor);
            iTwin = T2.TS(ir);
            
            iS = find(arrayfun(@(x) x.gID == ID_current,struCell{iE}));
            
            ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
            if isempty(ind)
                tGbVol = 0;
            else
                tGbVol = struCell{iE}(iS).tGbVol{iTwin}(ind);
            end
            
            % add to new colume
            T2.tGbVol(ir) = tGbVol;
            T2.tGbVolPct(ir) = tGbVol/(pi/4*(T2.gDia(ir)/umPerDp)^2);
            T2.tGbStrength(ir) = T2.tGbVolPct(ir)/T2.gb_length(ir);
        end
        
        save(['temp_results/twin_gb_summary_',num2str(iE),'.mat'],'T','T2','-append');
    end
end
%% next, need to go to the code for summary.















