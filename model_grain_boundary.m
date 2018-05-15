% function [gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID,x,y,stepSize)
% gb_dir{i} = 'horizontal' or 'vertical'
% gb_s_pt{i} = [id of points belonging to this grain boundary]
% pt_pos(i,:) = [x,y of this point]
% pt_s_gb{i} = [id of grain boundaries belongs to this point]
%
% chenzhe, 2018-05-14

function [gb_dir, gb_s_pt, pt_pos, pt_s_gb, tripleLookup] = model_grain_boundary(ID,x,y)

stepSize = y(2) - y(1);

[boundaryTF,~,neighborID,tripleTF,~] = find_one_boundary_from_ID_matrix(ID);
tripleLookup = [x(tripleTF>0),y(tripleTF>0)];

ind = boundaryTF>0;
gbPoints = [ID(ind), neighborID(ind), x(ind), y(ind)];

% sort the two grain IDs corresponding to gb points
t = gbPoints(:,[1,2]);
t = sort(t,2);
gbPoints(:,[1,2]) = t;

uniquePair = unique(gbPoints(:,[1,2]),'rows');    % unique grain pairs

% reset variables
gb_dir = [];
gb_s_pt = [];
pt_s_gb = [];
pt_pos = [];


% eliminate gb with less than 3 pts
longUniquePair = [];
for igb = 1:length(uniquePair)
    inds = (gbPoints(:,1)==uniquePair(igb,1))&(gbPoints(:,2)==uniquePair(igb,2));
    segPts = gbPoints(inds,[3,4]);     % points of this grain boundary segment
    if size(segPts,1) >= 3
        longUniquePair = [longUniquePair; uniquePair(igb,:)];
    end
end
uniquePair = longUniquePair;

ptCount = 1;    % count number of control points
for igb = 1:length(uniquePair)
    inds = (gbPoints(:,1)==uniquePair(igb,1))&(gbPoints(:,2)==uniquePair(igb,2));
    segPts = gbPoints(inds,[3,4]);     % points of this grain boundary segment
    
    % [1] determine end points of grain boundary. determine and record direction
    if range(segPts(:,1))>=range(segPts(:,2))
        segPts = sortrows(segPts,1);
        gb_dir{igb} = 'horizontal';
    else
        segPts = sortrows(segPts,2);
        gb_dir{igb} = 'vertical';
    end
    
    % determine up to 5 keypoints of this grain boundary segment for fitting
    inds = round(linspace(1,size(segPts,1), min(5, size(segPts,1))));
    keyPts = segPts(inds,:);
    
    for jj=1:size(keyPts,1)
        pt = keyPts(jj,:);  % current point considered
        
        % if first or last of keypoints, it could be a triple point
        if (jj==1)||(jj==size(keyPts,1))
            ind_triple = find(pdist2(pt,tripleLookup) <= 1.4143*stepSize,1,'first');
            if ~isempty(ind_triple)
                pt = tripleLookup(ind_triple,:);
            end
        end
        
        % determine if the point was (e.g., a triple point that was) already used before
        try
            [~,loc] = ismember(pt,pt_pos,'rows');
        catch
            loc = 0;
        end
        
        if loc > 0
            ipt = loc;
            ptCount = ptCount - 1;  % to conteract the ++ at the end of loop
        else
            ipt = ptCount;
        end
        
        % [2] record this grain boundary's point id
        try
            gb_s_pt{igb} = [gb_s_pt{igb}, ipt];
        catch
            gb_s_pt{igb} = ipt;
        end
        
        % [3] record this 'pt'
        pt_pos(ipt,:) = pt;
        
        % [4] record this point's grain boundary id
        try
            pt_s_gb{ipt} = [pt_s_gb{ipt}, igb];
        catch
            pt_s_gb{ipt} = igb;
        end
        
        % [at the end of loop] increment ipt
        ptCount = ptCount + 1;
        
    end
end

end

