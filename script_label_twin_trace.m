

% script to run a single iteration of function label_twin_trace
% function [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace(...
%     twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
%     struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,twinTF,debugTF,th_1,th_2, ssAllowed, goBack)

% Gather necessary info for script mode
goBack = 0;
th_1 = 0.3;
th_2 = 0.3;



if ~exist('goBack','var')
    goBack = 1;
end

% first check if the active system is the same. If so, skip the code
alreadyActive = struCell{iE_list(iEC)}(iS).cActiveSS(iC_list(iEC),:);
%%
if sum(alreadyActive(:)-ssAllowed(:)) ~= 0
    
    [~, ~, nss, ntwin, ~] = define_SS(sampleMaterial,'twin');
    
    % th_1; % max(traceVote) >= th_1 * length(peakAngles) to have enough_votes
    % th_2; % any traceVote > th_2 * max(traceVote) to become acceptable trace
    
    clusterNumMapL = clusterNumMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    if debugTF==1
        clusterNumMapL_a = clusterNumMapL;
        clusterNumMapL_a(ID_local~=ID_current)=nan;
        clusterNumMapL_a(clusterNumMapL_a==0)=nan;
        [f,a,c] = myplot(clusterNumMapL_a, grow_boundary(boundaryTFLocal));
        set(gca,'xticklabel','','yticklabel','');
        caxis([0 2]);
        set(c,'Ticks',[0 1 2]);
        title('Cluster ID','fontweight','normal');
        set(gca,'xTick',[],'yTick',[],'fontsize',18);
        colormap(parula(3));
        set(c,'limits',[2/3,2],'Ticks',[1,1+2/3],'TickLabels',{'1','2'});
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
    if debugTF==1
        myplot(clusterNumMapC, grow_boundary(boundaryTFLocal));
        title('cluster considered');
        set(gca,'xticklabel','','yticklabel','');
        title('');
    end
    
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
        % clusterNumMapT = double(thin(clusterNumMapC,inf));
        % clusterNumMapT = double(bwmorph((clusterNumMapC),'thin',inf) );
    end
    if debugTF==1
        % [f,a,c] = myplot(clusterNumMapT,  grow_boundary(boundaryTFLocal));
        t = clusterNumMapT;
        t(t==0) = 0.01;
        myplotm(t,boundaryTFLocal)
        ca([0, 2]);
        title('skeleton of cluster considered');
        set(gca,'xticklabel','','yticklabel','','fontsize',18);
        title('Skeleton', 'fontweight', 'normal');
        colorbar('off');
        
        myplotm(clusterNumMapL, 'TF',clusterNumMapT, 'r', 1);
        caxis([-0.1, max(clusterNumMapL(:))+0.1]);
        title('skeleton overlaid on cluster number map');
        set(gca,'xticklabel','','yticklabel','');
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
        title('Hough Transform of Skeleton','fontweight','normal')
        axis normal;
        hold on;
        for k = 1:size(peaks,1)
            xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
            plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
        end
        xlabel('\theta, degrees'); ylabel('\rho');
        set(gca,'fontsize',16,'xTick',[-90:45:90]);
    end
    
    % (8) Find lines. (A) Should keep gap small to prevent joining too many irrelavent parts. (2) MinLength should be decent
    % This is mainly to [illustrate] approximately, what are the peaks that were found ?
    if debugTF==1
        Lines = houghlines(clusterNumMapC, Theta, Rho, peaks, 'FillGap',5, 'MinLength', 0.1*min(size(clusterNumMapL)));
        myplotm(clusterNumMapL,'TF',grow_boundary(boundaryTFLocal), 'r', 1);
        caxis([-0.1, max(clusterNumMapL(:))+0.1]);
        hold on;
        for k = 1:length(Lines)
            xy = [Lines(k).point1; Lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','r');
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'.','LineWidth',2,'Color','r');
            plot(xy(2,1),xy(2,2),'.','LineWidth',2,'Color','r');
        end
        set(gca,'fontsize',18,'xTick',[],'yTick',[]);
        title(''); colorbar off;
        
        myplot(clusterNumMapT, grow_boundary(boundaryTFLocal));
        hold on;
        for k = 1:length(Lines)
            xy = [Lines(k).point1; Lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','r');
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'.','LineWidth',2,'Color','r');
            plot(xy(2,1),xy(2,2),'.','LineWidth',2,'Color','r');
        end
        set(gca,'fontsize',18,'xTick',[],'yTick',[]);
        title(''); colorbar off;
        
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
            branchPoints = bwmorph(skl, 'branchpoints');
            
            inds = find(branchPoints==1);
            [suby,subx] = ind2sub(size(branchPoints),inds);
            
            % can do this again to ensure to break skeleton, but not sure if its in general good --------------- 
            heavyClean = 1;
            if heavyClean
                branchPoints = imdilate(branchPoints, ones(3));
            end
            branch = skl;
            branch(branchPoints>0) = 0;
            
            % (11.2) assign an ID to each skeleton branch
            branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
            branchNumbered(~branch) = 0;
            branchNumbered(branchPoints>0) = 0;     % --- ? modify ??? 
            % and get the unique ID of the branches
            uniqueBranchNum = unique(branchNumbered(:));
            uniqueBranchNum(uniqueBranchNum==0)=[];
            % [illustrate] skeleton branches
            if debugTF==1
                % myplotm(mod(branchNumbered,5)+logical(branchNumbered));
                myplot(skl,grow_boundary(boundaryTFLocal));
                hold on;
                plot(subx,suby,'.g','markersize',24);
            end
            set(gca,'xTick',[],'yTick',[],'fontsize',18);
            title('');
            colorbar off;
            
            % (11.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.
            % Assign the ts/ss ID to the branches, which can be considered as grouped.
            branchGrouped = zeros(size(branchNumbered));        
            [~,fragments_N] = city_block(branchNumbered);
            for ib = 1:length(uniqueBranchNum)
                model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(traceND - branchND);
                dAngle(~activeSS) = inf;
                [~,ind] = min(dAngle);
                %  (*) 2019-01-22. Here we might have something to improve: if doesn't match well, maybe just assign zero? ----------------------------------------------------------  
                if 1&&(min(dAngle)>15)
                    ind = nan;
                end
                
                branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;   
                fragments_N(fragments_N==uniqueBranchNum(ib)) = uniqueBranchNum(ib) + ind*0;
            end
            
            % plot branchGroupped
            if debugTF >= 1
                branchGrouped_a = branchGrouped;
                branchGrouped_a(isnan(branchGrouped_a)) = 0;
                [f,a,c] = myplot(branchGrouped_a, grow_boundary(boundaryTFLocal)); 
                hold on;
                plot(subx,suby,'.g','markersize',24);
            end
            cluster_ID_colorbar(19,22,2);
            set(gca,'xTick',[],'yTick',[],'fontsize',18);
            set(c,'ticks',[19 22],'ticklabels',{'1','4'});
            title('');


            
            % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
            [~,fragmentsGrouped] = city_block(branchGrouped);
            fragmentsGrouped(clusterNumMapC==0) = 0;
            
            % (*) Use together with previous (*)
            fragmentsGrouped(isnan(fragmentsGrouped)) = 0;
            
            % [illustrate] the fragments
            if debugTF >= 1
                myplot(fragmentsGrouped, imdilate(branch,[0 0 0; 0 1 1; 0 0 0])); caxis([18,24]);
                [f,a,c] = myplot(fragmentsGrouped, grow_boundary(boundaryTFLocal)); caxis([18,24]);
            end
            cluster_ID_colorbar(19,22,2);
            set(gca,'xTick',[],'yTick',[],'fontsize',18);
            set(c,'ticks',[19 22],'ticklabels',{'1','4'});
            title('');
            
    end
    
    if sum(activeSS(:))
        haveActiveSS = 1;
    else
        haveActiveSS = 0;
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
    
    
    %% plot for one fragment: the hough line and gb intersection

%     close all;
    if 1
    
        %         [~,fragments_N] = city_block(branchNumbered);
        trueTwinMapL = evalin('base','trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);');
        uniqueBoundary_local = evalin('base','uniqueBoundary(indR_min:indR_max, indC_min:indC_max);');
        uniqueBoundary_local = imdilate(uniqueBoundary_local,ones(3));
        uniqueBoundary_local(ID_local~=ID_current)=0;
        
        trueTwinMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain

        % Can initialize field
        struCell{iE}(iS).tGb = cell(1,ntwin);
        struCell{iE}(iS).tGbPts = cell(1,ntwin);
        for ii = 1:ntwin
            struCell{iE}(iS).tGbPts{ii} = {[]};
        end
        dist_cr = 10;    % ------------------------------- This is a quite arbitrarily selected criterion  ---- -
        npts_cr = 5;    % # of points within a certain distance to a gb ponint, then this frag can be considered as touching gb
    
%         fragments_N(trueTwinMapL==0) = 0;   
%         fragments_N(ID_local~=ID_current) = 0;
        % Here get fragments_N from trueTwimMapL, by one_pass_label again
        fragments_N = one_pass_label_8(trueTwinMapL);   % [~,fragments_N] = city_block(branchNumbered);
        fragments_N(trueTwinMapL==0) = 0;
    
        
        fragNum = 17;    % fragNum = 17, is a good example 54545454545454545454545454545454545454545454545454545454545454545454545454 
        myplot(fragments_N==fragNum,boundaryTFLocal); hold on;
        bn = fragNum;
        fragMap = fragments_N;          % choice-1: branch_N
        fragMap(fragments_N~=bn) = 0;   % choice-1: branch_N
        % check length use x_local and y_local
        xrange = range(x_local(fragMap>0));
        yrange = range(y_local(fragMap>0));
        span = sqrt(xrange^2+yrange^2);
        spanEnough = span > min(size(y_local,1),size(y_local,2))*(y_local(2)-y_local(1))*0.10    % min(height,width) * scale * pct
        
        if true || spanEnough
            tsNum = mode(trueTwinMapL(fragMap>0));    % most frequent element, is the ts number of this branch  % choice-1: branch_N
            
            iTwin = tsNum - nss;
            
            theta_target = traceND(iTwin);
            [H,Theta,Rho] = hough(fragMap,'RhoResolution',1);
            % block other peaks, if any
            H(:,abs(Theta-theta_target)>2) = 0;
            peaks = houghpeaks(H, 1);
            peakRho = Rho(peaks(1));
            peakTheta = Theta(peaks(2));
            
            lines = houghlines(uniqueBoundary_local, Theta, Rho, peaks, 'FillGap',999999);   % Note that the lines can be empty
            
            % This calculates the trace dir
            frag_skl = double( bwskel(imbinarize(fragMap)));
            model = fitlm(x_local(frag_skl==1), y_local(frag_skl==1));
            branchND = atand(-1/model.Coefficients.Estimate(2));
            dAngle = abs(theta_target  - branchND);
            
            if (~isempty(lines) && (dAngle<15))
                xy = [lines(1).point1; lines(1).point2];

                % method-2: for each end point of houghline, calculate a distMap, and check if there is enough twin points within certain distance
                currentTwinMap = (trueTwinMapL==tsNum); % map of the twin system considered
                gd = sqrt(4*struCell{iE}(iS).gVol/pi);  % calculate the equivalent grain diameter, in # of data points
                
                for jj = 1:2
                    distMap = zeros(size(fragMap));
                    distMap(xy(jj,2),xy(jj,1)) = 1;
                    distMap = bwdist(distMap);
                    
%                     % criterion: if mean dist of twinned region within mantle to this intersection is in [0.3, 0.7] range
%                     mask = (currentTwinMap==1)&(distMap<gd*0.1);
%                     meanDist = nanmean(distMap(mask(:)==1));    % if it is a triangle/trapezoid, the mean height should always be 0.5, (similarly the distance to the base)
%                     if (meanDist>0.0333*gd)&&(meanDist<0.0666*gd)
                    % [alternatively] criterion for twin touch grain boundar: this fragment has > npts_cr ponints within distance dist_cr. 
                    dist_cr = 0.03 * gd;
                    npts_cr = 10;
                    if sum(sum( (distMap<dist_cr)&(fragMap>0) )) > npts_cr 
                        if debugTF
%                             mask_boundary = find_boundary_from_ID_matrix(distMap<gd*0.1);
%                             uniqueBoundary_local(mask_boundary==1) = 1;
%                             [f,a,c] = myplot(fragments_N==fragNum,uniqueBoundary_local); hold on;
                            
                            plot3(xy(:,1),xy(:,2),[100;100],'LineWidth',1,'Color','r');
                            % Plot beginnings and ends of lines
                            plot3(xy(jj,1),xy(jj,2),[100;100],'o','LineWidth',2,'Color','g');
                            title('');
                            set(gca,'xTick',[],'yTick',[],'fontsize',18);
                            colorbar off;
                        end
                        % This is a contact point. Record the position, and the twin system
                        gbNum = uniqueBoundary_local(xy(jj,2),xy(jj,1));  % get gbNum
                        ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                        
                        % convert xy to absolute coordinate using x_local, y_local
                        xcoord = x_local(xy(jj,2),xy(jj,1));
                        ycoord = y_local(xy(jj,2),xy(jj,1));
                        
                        if isempty(ind)
                            % append
                            struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                            iGb = length(struCell{iE}(iS).tGb{iTwin});
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                        else
                            % add
                            iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                        end
                    end
                end
                
            end
        end
    end
    
    %% plot fragments + gb intersection
    if 2 
        
        uniqueFragment_N = unique(fragments_N(:));
        uniqueFragment_N(uniqueFragment_N==0) = [];
        % for each twin skeleton, need to have a minimum length along the trace direction
        activeTS = sum(struCell{iE}(iS).cTrueTwin,1);
        
        dist_cr = 10;    % ------------------------------- This is a quite arbitrarily selected criterion  ---- -
        npts_cr = 5;    % # of points within a certain distance to a gb ponint, then this frag can be considered as touching gb
        
        if debugTF
            [f,a,c] = myplot(fragmentsGrouped,boundaryTFLocal); hold on;
        end
        for ii = 1:length(uniqueFragment_N)
            bn = uniqueFragment_N(ii); % bn = 54 for the example
            fragMap = fragments_N;          % choice-1: branch_N
            fragMap(fragments_N~=bn) = 0;   % choice-1: branch_N
            % check length use x_local and y_local
            xrange = range(x_local(fragMap>0));
            yrange = range(y_local(fragMap>0));
            try
                span = sqrt(xrange^2+yrange^2);
            catch
                span = 0;
            end
            spanEnough = span > min(size(y_local,1),size(y_local,2))*(y_local(2)-y_local(1))*0.10;    % min(height,width) * scale * pct
            
            if spanEnough
                tsNum = mode(trueTwinMapL(fragMap>0));    % most frequent element, is the ts number of this branch  % choice-1: branch_N
                
                iTwin = tsNum - nss;
                
                theta_target = traceND(iTwin);
                [H,Theta,Rho] = hough(fragMap,'RhoResolution',1);
                % block other peaks, if any
                H(:,abs(Theta-theta_target)>2) = 0;
                peaks = houghpeaks(H, 1);
                peakRho = Rho(peaks(1));
                peakTheta = Theta(peaks(2));
                
                lines = houghlines(uniqueBoundary_local, Theta, Rho, peaks, 'FillGap',999999);   % Note that the lines can be empty
                
                % This calculates the trace dir
                frag_skl = double( bwskel(imbinarize(fragMap)));
                model = fitlm(x_local(frag_skl==1), y_local(frag_skl==1));
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(theta_target  - branchND);
                
                if (~isempty(lines) && (dAngle<15))
                    xy = [lines(1).point1; lines(1).point2];

                    % method-2: for each end point of houghline, calculate a distMap, and check if there is enough twin points within certain distance
                    currentTwinMap = (trueTwinMapL==tsNum); % map of the twin system considered
                    gd = sqrt(4*struCell{iE}(iS).gVol/pi);  % calculate the equivalent grain diameter, in # of data points

                    for jj = 1:2
                        distMap = zeros(size(fragMap));
                        distMap(xy(jj,2),xy(jj,1)) = 1;
                        distMap = bwdist(distMap);
                        
%                         % criterion: if mean dist of twinned region within mantle to this intersection is in [0.3, 0.7] range
%                         mask = (currentTwinMap==1)&(distMap<gd*0.1);
%                         meanDist = nanmean(distMap(mask(:)==1));    % if it is a triangle/trapezoid, the mean height should always be 0.5, (similarly the distance to the base)
%                         if (meanDist>0.0333*gd)&&(meanDist<0.0666*gd)
                        % [alternatively] criterion for twin touch grain boundar: this fragment has > npts_cr ponints within distance dist_cr.
                        dist_cr = 0.03 * gd;
                        npts_cr = 10;
                        if sum(sum( (distMap<dist_cr)&(fragMap>0) )) > npts_cr
                    
                            if debugTF
%                                 plot3(xy(:,1),xy(:,2),[100;100],'LineWidth',1,'Color','r');     % This is the line
                                % Plot beginnings and ends of lines
                                plot3(xy(jj,1),xy(jj,2),100,'o','LineWidth',2,'Color','g');
%                                 text(xy(jj,1),xy(jj,2),num2str([xy(jj,1),xy(jj,2),bn]));
                            end
                            % This is a contact point. Record the position, and the twin system
                            gbNum = uniqueBoundary_local(xy(jj,2),xy(jj,1));  % get gbNum
                            ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                            
                            % convert xy to absolute coordinate using x_local, y_local
                            xcoord = x_local(xy(jj,2),xy(jj,1));
                            ycoord = y_local(xy(jj,2),xy(jj,1));
                            
                            if isempty(ind)
                                % append
                                struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                                iGb = length(struCell{iE}(iS).tGb{iTwin});
                                struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                            else
                                % add
                                iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                                struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
    title('');
    cluster_ID_colorbar(19,22,2);
    set(gca,'xTick',[],'yTick',[],'fontsize',18);
    set(c,'ticks',[19,22],'ticklabels',{'1','4'});
    
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

% end