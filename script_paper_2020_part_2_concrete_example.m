% For selected grain
% label theoretical trace direction (highlight active?)
% find the active twin system (by looking at theoretical trace direction)  
% find the intersecting neighbor/all neighbors
% find available/all deformation modes in neighbor
% find the mPrime_matrix

clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_by_Matlab_after_realign','choose a path of the saved processed data, or WS, or etc.'),'\'];
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
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gDiameter','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gDiameter','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
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

[newVariantFile, newVariantFilePath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map.mat','select the new result of dividing twin into variants');
load(fullfile(newVariantFilePath,newVariantFile),'struCell','trueTwinMapCell');

%% Find triple points
[boundary,~, ~, tripleTF, ~, indTriple, tripleIDs] = find_one_boundary_from_ID_matrix(ID);
xTriple = X(indTriple);
yTriple = Y(indTriple);
%%
iE = 3;
% close all;
data = [];

plotTF = 1;
pauseTF = 0;
umPerX = 360/4096;
umPerDp = 360/4096*5;
dToTriple_th = 5;       % eliminate intersection whose distance to triple point is smaller than this value
dToTriple_th_to_label = dToTriple_th;  % label if distance of intersection to triple point is smaller than this value

% Initialize table
variableNames = {'iE','ID','gDia','ID_neighbor','gDia_neighbor','TS','TSF','gb_length',...
    'incoming','iiE_each_twin','iiE_each_twin_at_this_boundary','intersection_to_triple','iiE_twins_at_this_boundary_nb',...
    'mPrime','rank_mPrime','ssn_nb','SF_nb',...
    'resB','rank_resB','ssn_nb_r','SF_nb_r',...
    'mPrime_wrtB','rank_mPrime_wrtB','ssn_nb_wrtB','SF_nb_wrtB', 'resB_wrtB','rank_resB_wrtB','ssn_nb_r_wrtB','SF_nb_r_wrtB',...
    'mPrime_wrtT','rank_mPrime_wrtT','ssn_nb_wrtT','SF_nb_wrtT', 'resB_wrtT','rank_resB_wrtT','ssn_nb_r_wrtT','SF_nb_r_wrtT',...
    'initiating', 'eMean_1','eMean_2','eMean_1_nb','eMean_2_nb', 'max_basal_SF','max_twin_SF','max_basal_SF_nb','max_twin_SF_nb',...
    'exz_ba','exz_pr','exz_py','exz_pyII','exz_etw', 'exzr_ba','exzr_pr','exzr_py','exzr_pyII','exzr_etw',...
    'tGbVol','tGbVolPct','tGbStrength'};
T_template = cell2table(cell(0,length(variableNames)));
T_template.Properties.VariableNames = variableNames;

% For a selected grain of interest
iS = find(arrayfun(@(x) x.gID == 425,struCell{iE})); 
ID_current = struCell{iE}(iS).gID;
ind = find(gID==ID_current);
euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
% g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
[abs_schmid_factor, ~, ~] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);

nNeighbors = gNNeighbors(ind);
ID_neighbors = gNeighbors(ind, 1:nNeighbors);

ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);

% Make it one data point wider on each side
indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1) + 100;
indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1) -300;
indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1) + 100;
indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1) -300;

ID_local = ID(indR_min:indR_max, indC_min:indC_max);
X_local = X(indR_min:indR_max, indC_min:indC_max);
Y_local = Y(indR_min:indR_max, indC_min:indC_max);
uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
exx_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);

% find out grain size, in um
gDia = sqrt(4*(struCell{iE}(iS).gVol* (umPerDp)^2)/pi);
% Find active system, if any, using cTrueTwin/tGb field
% activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
activeTS = cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb);
% convert to double with nan replacing 0
activeTS = double(activeTS);
activeTS(activeTS==0)=nan;
        
tSF = struCell{iE}(iS).tSF(:);

% If grain of interest is twinned (here it means: got a twin-gb intersection labeled)
if any(cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb)) % sum(struCell{iE}(iS).cTrueTwin(:))>0
    if plotTF==1
        [handleFig,aa,~] = myplot(X_local, Y_local, exx_local, imdilate(boundaryTF_local,[0 0 0; 1 1 0; 0 0 0]));  
%         [handleFig,aa,~] = myplotm(exx_local,'x',X_local, 'y',Y_local, 'tf', imdilate(boundaryTF_local,[0 0 0; 1 1 0; 0 0 0]));
        colormap(summer);
%         label_map_with_ID(X_local, Y_local, ID_local, handleFig, ID_current,'k',18);
        local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_current, [19,22], gca);
        disableDefaultInteractivity(aa);
        hold on;
    end
    
    % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
    for iNb = 1:nNeighbors
        ID_neighbor = ID_neighbors(iNb);
        iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
        
        ind = find(gID==ID_neighbor);
        euler_nb = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        
        if ~isempty(iS_neighbor)
            gDia_neighbor = sqrt(4*(struCell{iE}(iS_neighbor).gVol* (360/4096*5)^2)/pi);
            activeTS_nb = cellfun(@(x) ~isempty(x), struCell{iE}(iS_neighbor).tGb);
            % convert to double with nan replacing 0
            activeTS_nb = double(activeTS_nb);
            activeTS_nb(activeTS_nb==0)=nan;
                    
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
            
            [exz_iTwin_jMode,exz_iTwin_jSS] = calculate_exz(euler, euler_nb,sampleMaterial);
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

                % Based on whether considering basal or twin in this neighbor, do different things: ...
                % (3) label basal/twin slip trace
                ID_target = ID_neighbor;
                if plotTF==1
                    local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, outgoing_ssn_0, gca);
                    % special case --------------------------------------------------------------------------------------------------------------special----   
                    if ismember(ID_neighbor,386)
                        local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, [19:24], gca);
                    end
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
                    plot(X_local(inds),Y_local(inds),'.','markersize',18,'color',[1 0 1] * (mPrime_toPlot+1)/2);
%                     if outgoing_basal==1
%                         text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','r','fontsize',16);
%                     else
%                         text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','b','fontsize',16);
%                     end
                end
                
                % Extend table for summary.  -- can assign in local table using name?
                                
            else
                disp('not valid_grain_a and valid_grain_b')
            end
            % end of (valid_grain_a)&&(valid_grain_b)
            
        end
        % end of ~isempty(iS_neighbor)
        data(iNb).nb = ID_neighbor;
        data(iNb).mPrime = mPrimeMatrix(19:24,[1:6,19:24]);
        data(iNb).exz = exz_iTwin_jMode;
        data(iNb).exzss = exz_iTwin_jSS;
        data(iNb).SF = schmidFactorG1(19:24);
        data(iNb).SF_nb = schmidFactorG2([1:6,19:24]);
        data(iNb).euler = euler;
        data(iNb).euler_nb = euler_nb;
    end
    % end of for iNb = 1:nNeighbors
    



end

set(gca,'fontsize',16,'xTick',[],'yTick',[]);
title('');
% ii = 4;
% small_table = [];
% [data(ii).SF(:), data(ii).mPrime]

%% scripts useful
print('temp_results/AOI_iE=5.tiff','-r300','-dtiff')

%% Arrange text in one figure
a = findobj(gcf,'type','text');
%% clean in another figure, and copy  
b = findobj(gcf,'type','text');
delete(b);
copyobj(a,gca);
%%
hcp_cell('euler',data(2).euler_nb,'stress',[-1 0 0;0 0 0; 0 0 0],'ss',25:30)



