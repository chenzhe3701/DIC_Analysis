% function distMap = distance_from_boundary_in_grain(ID, gbList, gList, roundTF)
function distMap = distance_from_boundary_in_grain(ID, bgList, roundTF)
% Cauculate the (euclidean) distance map on the ID map, from
% unique_boundaries in gbList, in grains in gList
% Chenzhe, 2019-02-06

% this might be slow
% for debug:
% ID = [
%     1 1 1 1 1 2;
%     1 1 1 2 2 2;
%     1 1 1 2 2 2;
%     1 1 3 3 3 3;
%     1 1 3 3 3 3;
%     3 3 3 3 3 3;
%     ];
% gbList = [30002];
% gList = [3];

if ~exist('roundTF','var')
    % disp('distance round to integer');
    roundTF = 1;
end

% if length(gbList)~=length(gList)
%    error('gbList and gList must be of the same length'); 
% end

if ~isempty(bgList)
    gList = bgList(:,2);
    gbList = bgList(:,1);
else
    gList = [];
    gbList = [];
end
list = sortrows([gList(:),gbList(:)]);


[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);

distMap = zeros(size(ID));

method = 'slow';
switch method
    case 'slow'
        % slow method
        while ~isempty(list)
%         for ii = 1:length(gList)
            ID_current = list(1,1); % disp(ID_current);
            ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
            % Make it one data point wider on each side
            indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
            indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
            indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
            indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
            uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
            
            ind = (list(:,1)==ID_current);
            gbs_current = list(ind,2);
            list(ind,:) = [];   % remove the gbs of intereste related to this grain  
            
            BW_local = zeros(size(ID_local));   % initalize with 0
            BW_local(ismember(uniqueBoundary_local, gbs_current)) = 1;   % label initial boundary
            BW_local = bwdist(BW_local);    % calculate dist
            BW_local(~ismember(ID_local,ID_current)) = 0;
            
            distMap(indR_min:indR_max, indC_min:indC_max) = distMap(indR_min:indR_max, indC_min:indC_max) + BW_local;
        end
    otherwise
        [S, gID, gNNeighbors, gNeighbors] = generate_neighbor_structure_from_ID_map(ID);
        % algorithm: iteratively find grain that are not connected.  But there are still bugs, because a grain is still affected by non-connected close grains ...  
        while ~isempty(gbList)
            gList_cycle = [];
            gbList_cycle = [];
            
            for ii=1:length(gbList)
                gb = gbList(ii);
                g1 = floor(gb/10000);
                g2 = mod(gb,10000);
                
                if ismember(g1,gList) % if g1 is a grain of interest
                    g1_nbs = gNeighbors(gID==g1, 1:gNNeighbors(gID==g1));
                    g1_nbs(g1_nbs==g2) = [];
                    if ~any(ismember(g1_nbs,gList_cycle)) % if g1 does not already have a neighbor that is in the todolist in this cycle
                        gList_cycle = [gList_cycle;g1];
                        gbList_cycle = [gbList_cycle;gb];
                    end
                end
                
                if ismember(g2,gList) % if g2 is a grain of interest
                    g2_nbs = gNeighbors(gID==g2, 1:gNNeighbors(gID==g2));
                    g2_nbs(g1_nbs==g1) = [];
                    if ~any(ismember(g2_nbs,gList_cycle)) % if g1 does not already have a neighbor that is in the todolist in this cycle
                        gList_cycle = [gList_cycle;g2];
                        gbList_cycle = [gbList_cycle;gb];
                    end
                end
            end
            
            gList_cycle = unique(gList_cycle);
            gbList_cycle = unique(gbList_cycle);
            gbList(ismember(gbList,gbList_cycle)) = []; % remove grain boundary processed in this cycle from the list
            
            distMap_cycle = zeros(size(ID));
            distMap_cycle(ismember(uniqueBoundary, gbList_cycle)) = 1;
            valid_mask = ismember(ID,gList_cycle);
            distMap_cycle = bwdist(distMap_cycle);
            distMap_cycle(~valid_mask) = 0;
            
            distMap = distMap + distMap_cycle;
        end
end

if roundTF==1
    distMap = round(distMap);
end

end