
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

function [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace(...
    twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
    struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,twinTF,debugTF,th_1,th_2, ssAllowed, goBack)

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
    ok_1 = (pdistCS(:)<1.5*minPdist);
    
    eCluster = effective_strain_nBy3(struCell{iE}(iS).cCen(iC,:));
    eTwin = effective_strain_nBy3(struCell{iE}(iS).tStrain);
    ok_2 = (eCluster > 0.4*eTwin(:)) & (eCluster < 1.5*eTwin(:));
    
    cAllCluster = effective_strain_nBy3(struCell{iE}(iS).cCen);
    [sorted,rank_in_raw] = sort(cAllCluster);
    [tf,rank_in_sorted] = ismember(cAllCluster,sorted);
    rank_in_sorted_0_base = rank_in_sorted - 1; % convert to 0-based rank for easier comparison
    ok_3 = rank_in_sorted_0_base(iC) >= 0.3 * (length(struCell{iE}(iS).cLabel)-1);
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
    angleToBlock = 9;    % if use tracking, block within an angle range of the already active SS.
    % [To prevent it is an all-zero map].
    if sum(peakStrength)>0
        for ip = 1:length(peakAngles)
            
            % Can preset a threshold, but maybe just record the SF rather than using a threshold 
            SF_th = -0.5; 
            dAngle = abs(traceND - peakAngles(ip));
            dAngle(dAngle > angleThreshold) = inf;
            dAngle(dAngle < 1) = 1;
            % here we want to achieve that, for dAngle sasitfied, even if traceSF < 0, it still contributes
            traceSF_logsig = logsig(transfer_to_logsig(traceSF, 0.25, 0.4, 0.8));
            score = traceSF_logsig./dAngle;
            score(traceSF < SF_th) = 0;
            
            % [add someting] if there was already an activeSS, any trace within +-XXX degree will have reduced voting power, reduce score to 0.
            for ii=1:length(active_pre)
                if (active_pre(ii)==1)%&& (sum(abs(traceND-traceND(ii))<10) > 1)
                    ind = (abs(traceND-traceND(ii))>0)&(abs(traceND-traceND(ii))<angleToBlock);
                    score(ind) = 0;
                end
            end
            
            % normalize.  Need to check if max value > 0, becuase it can be all zero. 
            if max(score)>0
                score = score/max(score);
            end
            traceVote = traceVote + score;
        end
    end
    
    % [Need enough distinct peaks] The voted trace should be distinct. So if max(traceVote) < th_1 (e.g., 30%) * length(peakAngles), that means it's 'junk' vote
    enough_votes = max(traceVote) >= length(peakAngles) * th_1;
    
    traceOKSS = (traceVote > th_2 * max(traceVote)) .* enough_votes;    % Any one larger than 30% max vote is also selected --> this need re-tunning
    
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
            fragmentsGrouped = zeros(size(clusterNumMapC));
            %             fragmentsR2 = zeros(size(clusterNumMapC));
        case 1     % If try always run fit r2 procedure, can disable this part by using case '-1'
            ind = find(activeSS);
            fragmentsGrouped = (nss+ind) * ones(size(clusterNumMapC));
            fragmentsGrouped(clusterNumMapC==0) = 0;
            if debugTF >= 1
                myplot(fragmentsGrouped); caxis([18,24]);
            end
            %             fragmentsR2 = ones(size(clusterNumMapC));   % no other choices,  simply make it as one.
            %             fragmentsR2(clusterNumMapC==0) = 0;
        otherwise
            % (11) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
            % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
            
            % (11.1) break skeleton into small branches
            skl = clusterNumMapT;
            % branchPoints = bwmorph(skl, 'branchpoints');    % blocked, [chenzhe 2019-02-03]
            
            % can do this again to ensure to break skeleton, but not sure if its in general good ---------------
            heavyClean = 1;
            if heavyClean
                branchPoints = imdilate(branchPoints, ones(3));
            end
            branch = skl - branchPoints;    
            branch(branch<0)=0; % added, should it? [chenzhe 2019-02-03] 
            
            % (11.2) assign an ID to each skeleton branch
            branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
            branchNumbered(~branch) = 0;
            % and get the unique ID of the branches
            uniqueBranchNum = unique(branchNumbered(:));
            uniqueBranchNum(uniqueBranchNum==0)=[];
            % [illustrate] skeleton branches
            if debugTF==1
                myplotm(mod(branchNumbered,5)+logical(branchNumbered));
            end
            
            % (11.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.
            % Assign the ts/ss ID to the branches, which can be considered as grouped.
            branchGrouped = zeros(size(branchNumbered));
            for ib = 1:length(uniqueBranchNum)
                model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(traceND - branchND);
                for ii=1:length(dAngle)
                   if dAngle(ii) > 90
                       dAngle(ii) = 180 - dAngle(ii);   % this is how you findout the diff between angles ----------------------------------------  
                   end
                end                
                dAngle(~activeSS) = inf;
                [~,ind] = min(dAngle);
                %  (*) 2019-01-22. Here we might have something to improve: if doesn't match well, maybe just assign zero? ----------------------------------------------------------
                %  if set to >90, means just accept it, use closest match
                if 1&&(min(dAngle)>90)
                    ind = nan;
                end
                
                branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;
            end
            
            % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
            [~,fragmentsGrouped] = city_block(branchGrouped);
            fragmentsGrouped(clusterNumMapC==0) = 0;
            
            % (*) Use together with previous (*)
            fragmentsGrouped(isnan(fragmentsGrouped)) = 0;
            
            % [illustrate] the fragments
            if debugTF >= 1
                myplot(fragmentsGrouped, branch); caxis([18,24]);
            end
            
    end
    
    if sum(activeSS(:))
        haveActiveSS = 1;
    else
        haveActiveSS = 0;
    end
    disp(ssAllowed(:)');
    disp(active_post(:)');
    if (sum(ssAllowed(:)-active_post(:))~=0)&&(sum(~ssAllowed(:))+sum(active_post(:))~=0)
       error('error'); 
    end
    % check if need to go back.  If have activeSS, and not the first in the iE_list, go back.  If go back, only activeSS is allowed.
    % If go back, it may affect the previous strain levels maps.
    goBack = goBack;
    if (goBack)&&(haveActiveSS)&&(iEC~=1)
        display([' ---------------> go back: ' num2str(iE_list(iEC-1))]);
        [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
            struCell,iS,iE_list(iEC-1),iC_list(iEC-1),iE_list,iC_list,iEC-1,iE_stop,traceND,traceSF,sampleMaterial,twinTF,debugTF, 0.25, 0.7, activeSS);
    else
        haveActiveSS = 0;   % why set it to zero?
    end
    
else
    haveActiveSS = 0;
    fragmentsGrouped = [];
end

% Only update iE level.  Other levels were updated iteratively.
if ~isempty(fragmentsGrouped)
    twinMapCell_cluster{iE,iC} = fragmentsGrouped;
    
    sfMap = zeros(size(fragmentsGrouped));
    for it = 1:ntwin
        sfMap(fragmentsGrouped==it+nss) = traceSF(it);
    end
    sfMapCell_cluster{iE,iC} = sfMap;
end

end