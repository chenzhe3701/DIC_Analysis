% script to count typical number of peaks in the hough transform

peakInfo = [];
for iS = 1:length(struCell{2})
    ID_current = stru(iS).gID
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    for iE = iE_start:iE_stop
        
        clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
        clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
        
        nC = length(struCell{iE}(iS).cLabel);
        
        for iC = 1:nC
            cNum = struCell{iE}(iS).cLabel(iC);
            
            clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
            clusterNumMapC(clusterNumMapC~=iC) = 0;
            
            clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
            [clusterNumMapT, branchPoints] = clean_skl(clusterNumMapT, round(min(size(clusterNumMapC))*0.05));
            
            % (6) Then do hough transform. H = intensity on [rhos, thetas] map
            [H,Theta,Rho] = hough(clusterNumMapT,'RhoResolution',1);    % H: d-1 is rho, d-2 is theta.
            % this is just to examine typical number of peaks if not limit to 32
            if sum(H(:))>0
                npeaks = size(houghpeaks(H, 1e10, 'Threshold', 0.3 * max(H(:)), 'NHoodSize',[round_odd(0.05*min(size(clusterNumMapL))),5] ),1);
            else
                npeaks = 0;
            end
            peakInfo = [peakInfo; iS, ID_current, iE, iC, npeaks];
        end
    end
end

