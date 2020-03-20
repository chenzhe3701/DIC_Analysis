% For selected grain
% label theoretical trace direction (highlight active?)
% find the active twin system (by looking at theoretical trace direction)  
% find the intersecting neighbor/all neighbors
% find available/all deformation modes in neighbor
% find the mPrime_matrix

iE = 5;
close all;
data = [];

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

% For selected grain, 
iS = find(arrayfun(@(x) x.gID == 425,struCell{iE})); 

ID_current = struCell{iE}(iS).gID
ind = find(gID==ID_current);

euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];

g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
[abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
traceDir = abs_schmid_factor(nss+1:nss+ntwin,3)

nNeighbors = gNNeighbors(ind);
ID_neighbors = gNeighbors(ind, 1:nNeighbors);

ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);

% Make it one data point wider on each side
indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1) + 150;
indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1) -300;
indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1) + 200;
indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1) -300;

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
tSF = struCell{iE}(iS).tSF(:)

% If grain of interest is twinned (here it means: got a twin-gb intersection labeled)
if any(cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb)) % sum(struCell{iE}(iS).cTrueTwin(:))>0
    if plotTF==1
%         [handleFig,aa,~] = myplot(X_local, Y_local, exx_local, boundaryTF_local);  % boundaryTF_local_2 = find_boundary_from_ID_matrix(trueTwinMapLocal>0)|(boundaryTF_local);
        [handleFig,aa,~] = myplotm(exx_local,'x',X_local, 'y',Y_local, 'tf', boundaryTF_local);
        label_map_with_ID(X_local, Y_local, ID_local, handleFig, ID_current,'k',18);
%         local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_current, 19:24, gca);
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
            % activeTS_nb = sum(struCell{iE}(iS_neighbor).cTrueTwin,1)>0;
            activeTS_nb = cellfun(@(x) ~isempty(x), struCell{iE}(iS_neighbor).tGb);
            
            % (1.1) Calculate this_uniqueGB number.
            if ID_current > ID_neighbor
                gb = ID_current * 10000 + ID_neighbor;
            else
                gb = ID_neighbor * 10000 + ID_current;
            end
            gb_length = sum(uniqueBoundary_local(:)==gb);
            
            % strain calculation in area of interest.
            distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
            ind = (distMap_local>0)&(distMap_local<=50);
            eMean_1 =  nanmean(eMap_local(ind));   % mean
            
            distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_neighbor]);
            ind = (distMap_local>0)&(distMap_local<=50);
            eMean_1_nb =  nanmean(eMap_local(ind));   % mean
            
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
                [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler, euler_nb, stressTensor, [1 0 0], sampleMaterial, 'twin');
                
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
                
                mPrime_local = mPrimeMatrixAbs(19:24, outgoing_ssn);
                [mPrime_each_twin, ind] = max(mPrime_local,[],2);     % value of interests: m'
                [sorted_mPrime_each_twin,~] = sort(mPrime_each_twin,'descend');
                [~,rank_mPrime_each_twin] = ismember(mPrime_each_twin, sorted_mPrime_each_twin);    % m'-rank
                ssn_nb_for_each_twin = outgoing_ssn(ind);   % ssn in neighbor, for each ss/ts in the grain of interest
                SF_nb_for_each_twin = SFs_nb(ind);          % the SF of that selected ssn in neighbor, for each ss/ts in the grain of interest
                
                resB_local = resBurgersMatrixAbs(19:24, outgoing_ssn);
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
                        [mPrime_toPlot, ind] = max(mPrime_each_twin.*incoming_TS_TF(:));
                        mPrime_rank_toPlot = rank_mPrime_each_twin(ind);
                    else
                        [mPrime_toPlot, ind] = max(mPrime_each_twin.*activeTS(:));
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
                    plot(X_local(inds),Y_local(inds),'.','markersize',18,'color',[1 0 0] * mPrime_toPlot);
                    if outgoing_basal==1
                        text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','r','fontsize',16);
                    else
                        text(mean(X_local(inds)), mean(Y_local(inds)), 50, [num2str(mPrime_toPlot,'%.3f'), ', (',num2str(mPrime_rank_toPlot,1),')'], 'color','b','fontsize',16);
                    end
                end
                
                % Extend table for summary.  -- can assign in local table using name?
                                
            else
                disp('not valid_grain_a and valid_grain_b')
            end
            % end of (valid_grain_a)&&(valid_grain_b)
            
        end
        % end of ~isempty(iS_neighbor)
        data(iNb).mPrime = mPrimeMatrixAbs(19:24,1:3);
        data(iNb).exz = exz_iTwin_jMode;
        data(iNb).SF = schmidFactorG1(19:24)
        data(iNb).SF_nb = schmidFactorG2(1:3)
    end
    % end of for iNb = 1:nNeighbors
    



end


ii = 4;
small_table = [];
[data(ii).SF(:), data(ii).mPrime, 


