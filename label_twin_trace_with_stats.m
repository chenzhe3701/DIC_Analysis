
% This function returns a cell that contains the cluster map labeled by twin system number.
% for each cluster, its cLabel [iC] might be different depending on the strain level [iE].  So, refer to: [iE_list, iC_list].
% The selected/input [iE] is the 'strain label' entry point for the anlaysis for the selected/input [iC] (the cLable at the entry point)
% It runs iteratively, and analyze all necessary strain levels, but not always all strain levels.
% In the output, each cell contains TS map of the cluster at the analyzed [iE]
% For iEs that is not anlayzed, the cell contains [].
%
% It can change the cluster's twin label/map in the current and previous iEs history.
% The label is stored in the tMap_iEC_entry_cell{iE,iC}, where iE, iC can be determined from iE_list, iC_list, iEC.
% This function is actually update the tMap_cell which is of max dimension (iE_stop, iC_max) 
%
% chenzhe, 2019-12-25
% label with different tunning parameters explicitly expressed. 
% Note that this does not include a final clean up.

function [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace_with_stats(...
    twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local,uniqueBoundary_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
    struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceDir,traceND,traceSF,sampleMaterial,twinTF,debugTF,th_1,p, ssAllowed, goBack)

diffStrain_cr = p.diffStrain_cr;    % 0.07;         % diff in strain need to < 0.07;
rEffStrain = p.rEffStrain;          %[0.1, 1.9];    % ratio of cluster effective strain and twin effective strain should be within this range   
strainRank_cr = p.strainRank_cr;    %1;             % rank need to >=1 (in 0-based) 
SF_th = p.SF_th;                    %0.2;           % twin SF need to be larger than this  (or 0.1)

pctTotalPeaks_bwd_cr = p.pctTotalPeaks_bwd_cr ;     %0.1; % when going backward, need at least these pct of total number of peaks to be considered as active trace  (or 0.05, 0.15)  
% note that this is 'th_1' input outside of this function  

pctTotalPeaks_fwd_cr = p.pctTotalPeaks_fwd_cr ;     %0.3; % when going forward, need at least these pct of total number of peaks to be considered as active trace  (or 0.35)
% This is 'th_1' as an input inside this function, for iterratively go back

% Note th_2 is no longer used.


if ~exist('goBack','var')
    goBack = 1;
end

% first check if the active system is the same. If so, skip the code
alreadyActive = struCell{iE_list(iEC)}(iS).cActiveSS(iC_list(iEC),:);
if sum(alreadyActive(:)-ssAllowed(:)) ~= 0
    
    [~, ~, nss, ntwin, ~] = define_SS(sampleMaterial,twinTF);
    
    % th_1; % max(traceVote) >= th_1 * length(peakAngles) to have enough_votes
    % th_2; % any traceVote > th_2 * max(traceVote) to become acceptable trace
    
    clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    if debugTF==1
        myplot(clusterNumMapL);
    end
    
    % (1) Decide if any is twin cluster?  Check if the strain is reasonable compared to theoretical twin strain.
    %
    % (1.1) Strain should be close, e.g., dist < 1.5 min_dist
    % (1.2) Effective strain should also be similar in magnitude. Otherwise, it is just closest, but not similar
    % (1.3) Maybe, select the clusters with higher strains, e.g. 1-out-of-2, 2-out-of-3, 2-out-of-4, 3-out-of-5, ...
    % (1.4, 1.5) Is it active in a previous/later step.  This is useful when go back to a previous step and check  
    cNum = struCell{iE}(iS).cLabel(iC);
    indClusterLocal = (clusterNumMapL==cNum);
    
    pdistCS = pdist2(struCell{iE}(iS).cCen(iC,:), struCell{iE}(iS).tStrain);
    minPdist = min(pdistCS);
    ok_1 = pdistCS(:)<diffStrain_cr; %(pdistCS(:)<1.5*minPdist);
    
    eCluster = effective_strain_nBy3(struCell{iE}(iS).cCen(iC,:));
    eTwin = effective_strain_nBy3(struCell{iE}(iS).tStrain);
    ok_2 = (eCluster > rEffStrain(1)*eTwin(:)) & (eCluster < rEffStrain(2)*eTwin(:));
    
    cAllCluster = effective_strain_nBy3(struCell{iE}(iS).cCen);
    [sorted,rank_in_raw] = sort(cAllCluster);
    [tf,rank_in_sorted] = ismember(cAllCluster,sorted);     % this can make sure equal get equal rank, e.g., [.1 .1 .2] get rank [1 1 3]
    rank_in_sorted_0_base = rank_in_sorted - 1; % convert to 0-based rank for easier comparison
    ok_3 = (rank_in_sorted(iC)-1) > strainRank_cr * (length(struCell{iE}(iS).cLabel)-1);  %0.3 * (length(struCell{iE}(iS).cLabel)-1);
    ok_3 = ones(size(ok_1)) * ok_3;
    
    % Known to be active in later step. (If not every ss is allowed, then it must be a check coming back from a later step)  
    % When coming back, we might be able to relax requirement for ok_1 and ok_2 
    if iEC < length(iE_list)
        active_post = reshape(struCell{iE_list(iEC+1)}(iS).cActiveSS(iC_list(iEC+1),:), ntwin, 1);
    else
        active_post = zeros(ntwin,1);
    end
    %strainOKSS =  ok_3 & ((ok_1 & ok_2)|(active_post));
    
    % Know to be active in a previous step   
    if iEC > 1
        active_pre = struCell{iE_list(iEC-1)}(iS).cActiveSS(iC_list(iEC-1),:);
        active_pre = active_pre(:);
    else
        active_pre = zeros(ntwin,1);
    end

    
    if debugTF >= 1
        disp(['eCluster = ', num2str(eCluster)]);
        disp(['all clusters strain = ',num2str(cAllCluster')]);
        disp(['higher strain cluster ?, ok_3 = ', num2str(ok_3(1))]);
        disp(table(pdistCS(:), double(ok_1(:)),eTwin(:), double(ok_2(:)),...
            'variableNames',{'pdistCS','ok_1','eTwin','ok_2'}));
    end

    
    clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
    clusterNumMapC(clusterNumMapC~=iC) = 0;
    
    % (2) Check if this cluster is too close to the grain boundary ---------------------------------------------------------------      
    if struCell{iE}(iS).cToGbDist(iC,end) < 50
        % This means OK, cluster not too close to gb
        clusterNumMapC(clusterNumMapC==iC) = 0;
        disp(['cluster too close to gb: ID=',num2str(struCell{iE}(iS).gID),',iE=',num2str(iE),',iC=',num2str(iC)]);
    end
 
    
    % (5) Then Do thinning/skeleton. The bwskel() function can perform some prunning at the same time.
    % Default no pruning. Because it uses 8-connectivity in bwskel_(), prunning sometimes makes analysis worse.
    thinTF = 1;    % This is mainly for debug, but looks like we should always consider thinning/skeleton.
    if thinTF
        clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
        [clusterNumMapT, branchPoints] = clean_skl(clusterNumMapT, round(min(size(clusterNumMapC))*0.05));
        % clusterNumMapT = double(thin(clusterNumMapC,inf));
        % clusterNumMapT = double(bwmorph((clusterNumMapC),'thin',inf) );
    end
    if debugTF==1
        myplotm(clusterNumMapL, 'TF',clusterNumMapT, 'r', 1);
        caxis([-0.1, max(clusterNumMapL(:))+0.1]);
    end
    
    % (6) Then do hough transform. H = intensity on [rhos, thetas] map
    [H,Theta,Rho] = hough(clusterNumMapT,'RhoResolution',1);    % H: d-1 is rho, d-2 is theta.
    
    % (7) Find peaks. Set a [neighborhood size], d_width = 5% map size, d_angle = 5 deg.
    maxNumPeaks = 32;
    peaks = houghpeaks(H, maxNumPeaks, 'Threshold', 0.3 * max(H(:)), 'NHoodSize',[round_odd(0.05*min(size(clusterNumMapL))),5] );
    peakAngles = Theta(peaks(:,2));
    peakStrength = H(sub2ind(size(H),peaks(:,1),peaks(:,2)));
    if debugTF >= 1
        disp( table(peakAngles(:),peakStrength(:),'VariableNames',{'PeakAngles','PeakStrength'}) );
    end
    
    % This is just to [illustrate] where the peak is in the hough space
    if debugTF==1
        myplot(Theta,Rho,H);
        xlabel('theta'); ylabel('rho');
        axis normal;
        hold on;
        for k = 1:size(peaks,1)
            xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
            plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
        end
    end
    
    % (8) Find lines. (A) Should keep gap small to prevent joining too many irrelavent parts. (2) MinLength should be decent
    % This is mainly to [illustrate] approximately, what are the peaks that were found ?
    if debugTF==1
        lines = houghlines(clusterNumMapC, Theta, Rho, peaks, 'FillGap',5, 'MinLength', 0.1*min(size(clusterNumMapL)));
        myplotm(clusterNumMapL, 'r', 1);
        caxis([-0.1, max(clusterNumMapL(:))+0.1]);
        hold on;
        for k = 1:length(lines)
            xy = [lines(k).point1; lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','g');
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','r');
            plot(xy(2,1),xy(2,2),'^','LineWidth',2,'Color','r');
        end
    end
    
    % (9) Determine the active variant/ss by matching the peakAngles with traceND.
    % Use a [5 deg] threshold. Then, for those within valid angle range, make a score = SF/deltaAngle.
    traceVote = zeros(size(traceND));
    angleThreshold = 5;
    % angleToBlock = 9;    % if use tracking, block within an angle range of the already active SS.
    % [To prevent it is an all-zero map].
    if sum(peakStrength)>0
        for ip = 1:length(peakAngles)
            % Can preset a threshold, but maybe just record the SF rather than using a threshold
            dAngle = min( abs(traceND - peakAngles(ip)), abs(traceND + 180*(-1)*sign(traceND)- peakAngles(ip)) );
            dAngle = ceil(dAngle);
            % if multiple traces within angle of 5,  select the one with higher SF  
            if sum(dAngle<=angleThreshold)>1
                ind = (dAngle<=angleThreshold);
                maxSF = max(traceSF(ind));
                ind = (traceSF==maxSF);
                dAngle(~ind) = dAngle(~ind) + angleThreshold;
            end
            
            for itwin = 1:6
                if dAngle(itwin)<=angleThreshold
                    traceVote(itwin) = traceVote(itwin) + 1;
                end
            end
        end
    end
    
    
    
    % [Need enough distinct peaks] The voted trace should be distinct. So if max(traceVote) < th_1 (e.g., 30%) * length(peakAngles), that means it's 'junk' vote
    % enough_votes = max(traceVote) >= length(peakAngles) * th_1;
    % traceOKSS = (traceVote > th_2 * max(traceVote)) .* enough_votes;    % Any one larger than 30% max vote is also selected --> this need re-tunning
    traceOKSS =  (traceVote >= length(peakAngles) * th_1) & (traceSF>SF_th);
    
    % something left for debugging ...
    cVolPct = sum(clusterNumMapL(:)==iC)/sum(clusterNumMapL(:));
    
    % ----------------------- criterion -----------------------------------------------------------------------------------------------------------------------------------------
    % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % traceMatch
    % ok_1: strain close enough
    % ok_2: strain high enough
    % ok_3: it is the high strain cluster
    % refActiveSS: previously strain level active
    % ok_4: previous strain NOT active
    % cVolPctNotDecrease: cluster vol not decrease compared to previous strain level
    %     activeSS = ssAllowed & (traceOKSS & (ok_3 & ((ok_1&ok_2)|ok_4))) | refActiveSS(:);
    activeSS = ssAllowed & (traceOKSS & ok_3 & ((ok_1&ok_2)|active_post)) | active_pre(:);
    
    if debugTF >= 1
        disp(['# of peaks found: ', num2str(length(peakAngles))]);
        disp(['cluster vol pct: ', num2str(cVolPct)]);
        %     disp(['cluster vol not decrease: ', num2str(cVolPctNotDecrease)]);
        disp(table(traceSF,traceND,traceVote,double(traceOKSS),double(ok_1),double(ok_2),double(ok_3),active_pre(:),active_post(:),double(activeSS),...
            'variableNames',{'traceSF','traceND','traceVote','traceOk','ok_1','ok_2','ok_3','activePre','activePost','activeSS'}));
    end
    
    % record the activeSS
    struCell{iE}(iS).cActiveSS(iC,:) = activeSS;
    
    
    % If there are more than one active slip/twin systems, we should seperate them:
    switch sum(activeSS)
        case 0
            % fragmentsGrouped = zeros(size(clusterNumMapC));
            tnMap = zeros(size(clusterNumMapC));
        case 1     % If try always run fit r2 procedure, can disable this part by using case '-1'
            variant_num = find(activeSS);
            tnMap = zeros(size(clusterNumMapC));
            tnMap(clusterNumMapC==iC) = variant_num;
            if debugTF >= 1
                myplot(tnMap); caxis([1,6]);
            end
        otherwise
            clear csl;
            % need to divide
            for iTwin = 1:6
                if activeSS(iTwin)==1
                    % (Step-2) Rotate maps
                    ID_r = imrotate(ID_local, traceDir(iTwin), 'nearest', 'loose');
                    X_r = imrotate(x_local, traceDir(iTwin), 'nearest', 'loose');
                    Y_r = imrotate(y_local, traceDir(iTwin), 'nearest', 'loose');
                    %                         uniqueBoundary_r = imrotate(uniqueBoundary_local, traceDir(iTwin), 'bilinear', 'loose');
                    uniqueBoundary_r = imrotate(uniqueBoundary_local, traceDir(iTwin), 'nearest', 'loose');
                    uniqueBoundary_r = imdilate(uniqueBoundary_r,ones(3));
                    clusterNumMap_r = imrotate(clusterNumMapC, traceDir(iTwin), 'nearest', 'loose');
                    vMap_r = (clusterNumMap_r == iC); % try this to segment clusterNumMap into variantMap(or trueTwinMap)
                    
                    [nr,nc] = size(vMap_r);
                    gbLabelMap = zeros(nr,nc);    % to store assigned gb_label
                    gbNumXY_intersect = [];     % [gbNum, Xpos, Ypos]  ----------------------------------> this could be recorded in struCell.
                    cslMap = zeros(nr,nc);     % a map recording Connected Segment Length (CSL)   ----------> this might be helpful do determine variant number, if starting point is clusterNumMap.
                    gbLR = zeros(nr,2);     % store each rows two possible gbs
                    
                    % (Step-3)
                    for ir = 1:nr
                        if any(vMap_r(ir,:))
                            icL_back = find(uniqueBoundary_r(ir,:),1,'first');
                            gbL = uniqueBoundary_r(ir,icL_back);
                            icR_back = find(uniqueBoundary_r(ir,:),1,'last');
                            gbR = uniqueBoundary_r(ir,icR_back);
                            gbLR(ir,:) = [gbL, gbR];
                            
                            % (instert Step-4) determine if gbL/R can be considered as an intersecting gb. -----------------------------------
                            length_cr = round(min(30, (icR_back - icL_back)/2));
                            num_cr = round(length_cr * 0.7);
                            if sum(vMap_r(ir,icL_back:icL_back+length_cr))>num_cr
                                gbNumXY_intersect = [gbNumXY_intersect; gbL, X_r(ir,icL_back), Y_r(ir,icL_back)];
                            end
                            if sum(vMap_r(ir,icR_back-length_cr:icR_back))>num_cr
                                gbNumXY_intersect = [gbNumXY_intersect; gbR, X_r(ir,icR_back), Y_r(ir,icR_back)];
                            end
                            % end of (Step-4). Alternatively, read from manual label for data analysis. ---------------------------------------
                            
                            icL_front = find(vMap_r(ir,:),1,'first');
                            icR_front = find(vMap_r(ir,:),1,'last');
                            % will be false, if either is empty
                            while (icL_front<=icR_front)
                                csl_length = 0;
                                if (icL_front-icL_back)<=(icR_back-icR_front)
                                    % search for connected segments from left to right
                                    while(vMap_r(ir,icL_front))
                                        gbLabelMap(ir,icL_front) = 1;  % prepare to assign label
                                        vMap_r(ir,icL_front) = 0; % make element on variant map 0
                                        icL_front = icL_front + 1;   % move pointer forward to the right
                                        csl_length = csl_length + 1;
                                    end
                                    cslMap(gbLabelMap==1) = csl_length;
                                    gbLabelMap(gbLabelMap==1) = gbL;
                                    icL_back = icL_front - 1;    % assign left side back
                                    icL_front = find(vMap_r(ir,:),1,'first'); % search left side front again
                                else
                                    while (vMap_r(ir,icR_front))
                                        gbLabelMap(ir,icR_front) = 1;
                                        vMap_r(ir,icR_front) = 0;
                                        icR_front = icR_front - 1;
                                        csl_length = csl_length + 1;
                                    end
                                    cslMap(gbLabelMap==1) = csl_length;
                                    gbLabelMap(gbLabelMap==1) = gbR;
                                    icR_back = icR_front + 1;
                                    icR_front = find(vMap_r(ir,:),1,'last');
                                end
                            end
                        end % end of if any(variant_r(ir,:))
                    end % end of for ir=1:nr
                    
                    % (Step-5) go back to clean.  Sometimes, no intersection is determined ...
                    cleanTF = 1;
                    if (cleanTF)&&(~isempty(gbNumXY_intersect))
                        gbList = unique(gbNumXY_intersect(:,1));
                        for ir=1:nr
                            tf = ismember(gbLR(ir,:),gbList);
                            if sum(tf)==1
                                gbOK = gbLR(ir, tf);
                                gbKO = gbLR(ir, ~tf);
                                ind = gbLabelMap(ir,:) == gbKO;
                                gbLabelMap(ir,ind) = gbOK;
                            elseif sum(ismember(gbLR(ir,:),gbList))==0
                                gbLabelMap(ir,:) = 0;
                            end
                        end
                    end
                    
                    temp = cslMap;
                    % rotate back, need to crop again.
                    temp = imrotate(cslMap,-traceDir(iTwin), 'nearest', 'loose');
                    ID_back = imrotate(ID_r,-traceDir(iTwin), 'nearest', 'loose');
                    ind_back = ismember(ID_back, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
                    
                    % Need to crop a region of the original size.
                    img1_template = (ID_local==ID_current);
                    img2_signal = (ID_back==ID_current);
                    [yOffSet, xOffSet] = normxcorr2A_register(img1_template, img2_signal, [0 0 0 0], [0 0 0 0], 0);
                    indC_back_min = 1 + xOffSet;
                    indC_back_max = indC_back_min + indC_max - indC_min;
                    indR_back_min = 1 + yOffSet;
                    indR_back_max = indR_back_min + indR_max - indR_min;
                    
                    csl(:,:,iTwin) = temp(indR_back_min:indR_back_max, indC_back_min:indC_back_max);
                    
                end % end of if(activeTS(iTwin)==1)
                
            end % end of for iTwin=1:6
            
            [nr,nc,np] = size(csl);
            tnMap = zeros(nr,nc);   % variant_num_map of this cluster
            for ir=1:nr
                for ic=1:nc
                    [maxV,temp] = max(csl(ir,ic,:));
                    if maxV>0
                        tnMap(ir,ic) = temp;
                    end
                end
            end

            % note: 'fragmentsGrouped' = 'tnMap' in new method.  
            % [illustrate] the fragments
            if debugTF >= 1
                myplot(tnMap); caxis([1,6]);
            end
            
    end
    
    if sum(activeSS(:))
        haveActiveSS = 1;
    else
        haveActiveSS = 0;
    end
%     disp(ssAllowed(:)');
%     disp(active_post(:)');
    if (sum(ssAllowed(:)-active_post(:))~=0)&&(sum(~ssAllowed(:))+sum(active_post(:))~=0)
       error('error'); 
    end
    % check if need to go back.  If have activeSS, and not the first in the iE_list, go back.  If go back, only activeSS is allowed.
    % If go back, it may affect the previous strain levels maps.
    goBack = goBack;
    if (goBack)&&(haveActiveSS)&&(iEC~=1)
        display([' ---------------> go back: ' num2str(iE_list(iEC-1))]);
        [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace_with_stats(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local,uniqueBoundary_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
            struCell,iS,iE_list(iEC-1),iC_list(iEC-1),iE_list,iC_list,iEC-1,iE_stop,traceDir,traceND,traceSF,sampleMaterial,twinTF,debugTF, pctTotalPeaks_bwd_cr, p, activeSS);
    else
        haveActiveSS = 0;   % why set it to zero?
    end
    
else
    haveActiveSS = 0;
    tnMap = [];  % note: 'fragmentsGrouped' = 'tnMap' in new method.
end

% Only update iE level.  Other levels were updated iteratively.
if ~isempty(tnMap)
    twinMapCell_cluster{iE,iC} = tnMap;
    
    sfMap = zeros(size(tnMap));
    for it = 1:ntwin
        sfMap(tnMap==it) = traceSF(it);     % depend on if twin numbered from 1:6 or from 19:24, use tnMap==it or tnMap==it+nss  
    end
    sfMapCell_cluster{iE,iC} = sfMap;
end

end