
function [fragments, struCell, haveActiveSS] = label_twin_trace(cluster_number_maps_cleaned,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
    struCell,iS,iE,iC,iE_list,iC_list,iEC,traceND,traceSF,sampleMaterial,twinTF,debugTF,th_1,th_2, ssAllowed)

% first check if the active system is the same. If so, skip the code
alreadyActive = struCell{iE_list(iEC)}(iS).cActiveSS(iC_list(iEC),:);
if sum(alreadyActive(:)-ssAllowed(:)) ~= 0
    
    [~, ~, nss, ntwin, ~] = define_SS(sampleMaterial,twinTF);
    
    % th_1; % max(traceVote) >= th_1 * length(peakAngles) to have enough_votes
    % th_2; % any traceVote > th_2 * max(traceVote) to become acceptable trace
    
    clusterNumMapL = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    if debugTF==1
        myplot(clusterNumMapL);
    end
    
    % (4) Decide if any is twin cluster?  Check if the strain is reasonable compared to theoretical twin strain.
    %
    % (4.1) Strain should be close, e.g., dist < 1.5 min_dist
    % (4.2) Effective strain should also be similar in magnitude. Otherwise, it is just closest, but not similar
    % (4.3) Maybe, select the clusters with higher strains, e.g. 1-out-of-2, 2-out-of-3, 2-out-of-4, 3-out-of-5, ...
    % (4.4) Not active before
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
    
    %strainOKSS =  ok_1 & ok_2 & ok_3;
    
    if debugTF >= 1
        disp(['eCluster = ', num2str(eCluster)]);
        disp(['all clusters strain = ',num2str(cAllCluster')]);
        disp(['higher strain cluster ?, ok_3 = ', num2str(ok_3(1))]);
        disp(table(pdistCS(:), double(ok_1(:)),eTwin(:), double(ok_2(:)),...
            'variableNames',{'pdistCS','ok_1','eTwin','ok_2'}));
    end
    clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
    clusterNumMapC(clusterNumMapC~=iC) = 0;
    cVolPct = sum(clusterNumMapC>0)/sum(clusterNumMapL>0);
    % if cVolPct < cVolPctOld
    %     cVolPctNotDecrease = 0;
    % else
    %     cVolPctNotDecrease = 1;
    % end
    % cVolPctOld = cVolPct;
    % myplot(clusterNumMapC);
    
    % (5) Then Do thinning/skeleton. The bwskel() function can perform some prunning at the same time.
    % Default no pruning. Because it uses 8-connectivity in bwskel_(), prunning sometimes makes analysis worse.
    thinTF = 1;    % This is mainly for debug, but looks like we should always consider thinning/skeleton.
    if thinTF
        clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
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
    
    % [10] Should also consider the activeSS from previous step, and combine
    if iEC > 1
        refActiveSS = struCell{iE_list(iEC-1)}(iS).cActiveSS(iC_list(iEC-1),:);
    else
        refActiveSS = zeros(ntwin,1);
    end
    ok_4 = ~refActiveSS(:);
    
    % (9) Determine the active variant/ss by matching the peakAngles with traceND.
    % Use a [5 deg] threshold. Then, for those within valid angle range, make a score = SF/deltaAngle.
    traceVote = zeros(size(traceND));
    angleThreshold = 5;
    % [To prevent it is an all-zero map].
    if sum(peakStrength)>0
        for ip = 1:length(peakAngles)
            dAngle = abs(traceND - peakAngles(ip));
            dAngle(dAngle > angleThreshold) = inf;
            dAngle(dAngle < 1) = 1;
            %             score = traceSF./dAngle;  % here we want to achieve that, for dAngle sasitfied, even if traceSF < 0, it still contributes
            score = logsig(transfer_to_logsig(traceSF, 0.2, 0.4, 0.9)) ./ dAngle;
            
            % [add someting] if there was already an activeSS, any trace within +-10 degree will have reduced voting power, reduce score to 0.
            for ii=1:length(refActiveSS)
                if (refActiveSS(ii)==1)%&& (sum(abs(traceND-traceND(ii))<10) > 1)
                    %                 for jj = 1:length(traceND)
                    %                     if (abs(traceND(jj)-traceND(ii))>0)&&(abs(traceND(jj)-traceND(ii))<10)
                    %                         score(jj) = score(jj)/2*0;
                    %                     end
                    %                 end
                    ind = (abs(traceND-traceND(ii))>0)&(abs(traceND-traceND(ii))<12);
                    score(ind) = 0;
                end
            end
            
            % normalize
            if max(score)>0
                score = score/max(score);
            end
            traceVote = traceVote + score;
        end
    end
    
    
    
    
    % [Need enough distinct peaks] The voted trace should be distinct. So if max(traceVote) < th_1 (e.g., 30%) * length(peakAngles), that means it's 'junk' vote
    enough_votes = max(traceVote) >= length(peakAngles) * th_1;
    
    traceOKSS = (traceVote > th_2 * max(traceVote)) .* enough_votes;    % Any one larger than 30% max vote is also selected --> this need re-tunning
    
    %     % [Additionally] If traces match super good, but clusterSize very small, then maybe it's ok. --------------------> This criterion need tunning.
    %     % Where is the cluster 'small' criterion? .
    cVolPct = sum(clusterNumMapL(:)==iC)/sum(clusterNumMapL(:));
    %     small_cluster_good_trace = zeros(size(traceOKSS));
    %     [val, ind] = max(traceVote);
    %     if val >= 0.75 * length(peakAngles)
    %         small_cluster_good_trace(ind) = 1;
    %     end
    
    %     activeSS = strainOKSS & traceOKSS | small_cluster_good_trace;  % combine strainOK and traceOK
    
    
    % ----------------------- criterion -----------------------------------------------------------------------------------------------------------------------------------------
    % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % traceMatch
    % ok_1: strain close enough
    % ok_2: strain high enough
    % ok_3: it is the high strain cluster
    % refActiveSS: previously strain level active
    % ok_4: previous strain NOT active
    % cVolPctNotDecrease: cluster vol not decrease compared to previous strain level
    activeSS = ssAllowed & (traceOKSS & (ok_3 & ((ok_1&ok_2)|ok_4))) | refActiveSS(:);
    
    if debugTF >= 1
        disp(['# of peaks found: ', num2str(length(peakAngles))]);
        disp(['cluster vol pct: ', num2str(cVolPct)]);
        %     disp(['cluster vol not decrease: ', num2str(cVolPctNotDecrease)]);
        disp(table(traceSF,traceND,traceVote,double(traceOKSS),double(ok_1),double(ok_2),double(ok_3),double(ok_4),refActiveSS(:),double(activeSS),...
            'variableNames',{'traceSF','traceND','traceVote','traceOk','ok_1','ok_2','ok_3','ok_4','refActiveSS','activeSS'}));
    end
    
    % record the activeSS
    struCell{iE}(iS).cActiveSS(iC,:) = activeSS;
    
    % If there are more than one active slip/twin systems, we should seperate them:
    switch sum(activeSS)
        case 0
            fragments = zeros(size(clusterNumMapC));
        case 1
            ind = find(activeSS);
            fragments = (nss+ind) * ones(size(clusterNumMapC));
            fragments(clusterNumMapC==0) = 0;
            if debugTF >= 1
                myplot(fragments); caxis([18,24]);
            end
        otherwise
            % (11) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
            % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
            
            % (11.1) break skeleton into small branches
            skl = clusterNumMapT;
            branchPoints = bwmorph(skl, 'branchpoints');
            
            % can do this again to ensure to break skeleton, but not sure if its in general good ---------------
            heavyClean = 1;
            if heavyClean
                branchPoints = imdilate(branchPoints, ones(3));
            end
            branch = skl - branchPoints;
            
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
            % store the twin system r2 fit
            for itwin = 1:ntwin
                tR2{itwin} = 0;
            end
            for ib = 1:length(uniqueBranchNum)
                
                model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(traceND - branchND);
                dAngle(~activeSS) = inf;
                [~,ind] = min(dAngle);
                branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;
                
                r2 = model.Rsquared.Ordinary;
                if isnan(r2)
                    r2 = 0;
                end
                tR2{ind} = [tR2{ind},r2];
            end
            for itwin = 1:ntwin
                struCell{iE}(iS).tR2(iC,itwin) = mean(tR2{itwin});
            end
            %         if debugTF
            %             struCell{iE}(iS).tR2(iC,:)
            %         end
            % tR2 only have number when multiple ts are fitted.
            
            
            % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
            [~,fragments] = city_block(branchGrouped);
            fragments(clusterNumMapC==0) = 0;
            % [illustrate] the fragments
            if debugTF >= 1
                myplot(fragments, branch); caxis([18,24]);
            end
            
    end
    
    if sum(activeSS(:))
        haveActiveSS = 1;
    else
        haveActiveSS = 0;
    end
    
    % check if need to go back.  If have activeSS, and not the first in the iE_list, go back.  If go back, only activeSS is allowed.
    if (haveActiveSS)&&(iEC~=1)
        display(iEC);
        [~, struCell, haveActiveSS] = label_twin_trace(cluster_number_maps_cleaned,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
            struCell,iS,iE_list(iEC-1),iC_list(iEC-1),iE_list,iC_list,iEC-1,traceND,traceSF,sampleMaterial,twinTF,debugTF, 0.25, 0.7, activeSS);
    else
        haveActiveSS = 0;
    end
    
    
else
    haveActiveSS = 0;
    fragments = [];
end
