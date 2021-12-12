
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

function [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = modify_twin_trace(...
    twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
    struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,twinTF,debugTF,th_1,th_2, ssAllowed)



    [~, ~, nss, ntwin, ~] = define_SS(sampleMaterial,twinTF);
        
    clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
    clusterNumMapC(clusterNumMapC~=iC) = 0;
    % If cluster not too close to gb
    if struCell{iE}(iS).cToGbDist(iC,end) < 50
        clusterNumMapC(clusterNumMapC==iC) = 0;  
    end
    % If thinTF==1
    clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
    
    % record the activeSS, modify both fields.
    struCell{iE}(iS).cActiveSS(iC,:) = activeSS;
    
    % ------------------------------------------------------> what about .cTrueTwin(), and tVol() ?          
    struCell{iE}(iS).cTrueTwin(iC,:) = activeSS;
    
    % If there are more than one active slip/twin systems, we should seperate them:
    switch sum(activeSS)
        case 0
            fragments = zeros(size(clusterNumMapC));
        
        case 1     % If try always run fit r2 procedure, can disable this part by using case '-1'
            ind = find(activeSS);
            fragments = (nss+ind) * ones(size(clusterNumMapC));
            fragments(clusterNumMapC==0) = 0;

        otherwise
            % (11) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
            % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
            
            % (11.1) break skeleton into small branches
            skl = clusterNumMapT;
            branchPoints = bwmorph(skl, 'branchpoints');
            
            % If heavyClean = 1;
            branchPoints = imdilate(branchPoints, ones(3));
            branch = skl - branchPoints;
            
            % (11.2) assign an ID to each skeleton branch
            branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
            branchNumbered(~branch) = 0;
            % and get the unique ID of the branches
            uniqueBranchNum = unique(branchNumbered(:));
            uniqueBranchNum(uniqueBranchNum==0)=[];
            % [illustrate] skeleton branches
            % myplotm(mod(branchNumbered,5)+logical(branchNumbered));
            
            % (11.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.
            % Assign the ts/ss ID to the branches, which can be considered as grouped.
            branchGrouped = zeros(size(branchNumbered));
            for ib = 1:length(uniqueBranchNum)
                model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(traceND - branchND);
                dAngle(~activeSS) = inf;
                [~,ind] = min(dAngle);
                branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;
            end
            
            % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
            [~,fragments] = city_block(branchGrouped);
            fragments(clusterNumMapC==0) = 0;

            % [illustrate] the fragments
            % myplot(fragments, branch); caxis([18,24]);
            
    end


    % Only update iE level.  Other levels were updated iteratively.
    if ~isempty(fragments)
        twinMapCell_cluster{iE,iC} = fragments;
        
        sfMap = zeros(size(fragments));
        for it = 1:ntwin
            sfMap(fragments==it+nss) = traceSF(it);
        end
        sfMapCell_cluster{iE,iC} = sfMap;
    end

end