function [S, gID, gNNeighbors, gNeighbors] = generate_neighbor_structure_from_ID_map(ID,min_valid_gID)
% find out the neighbors of each grain from an ID map
% [s=neighborStruct] construct grain neighbor structure S.g1 = [1;2;3], S.g2{ind_in_g1} = [2;3;...]
% chenzhe, 2019-02-08

if ~exist('min_valid_gID','var')
    disp('mininum valid grain ID number set to 2');
    min_valid_gID = 2;  % e.g., grain # 1 represents the frame of the EBSD area, so it is not a valid grain.
end

[boundaryTF, boundaryID, neighborID, ~, ~,] = find_one_boundary_from_ID_matrix(ID);
pairs = unique([ID(boundaryTF==1), neighborID(boundaryTF==1)], 'rows');
ind = pairs(:,1)<min_valid_gID | pairs(:,2)<min_valid_gID;
pairs(ind,:) = [];
pairs = unique([pairs; fliplr(pairs)], 'rows');

clear gNNeighbors_aligned gNeighbors_aligned neighborStruct_aligned;


gID = unique(ID(:));
for ii = 1:length(gID)
    S.g1(ii) = gID(ii);
    
    % find out how many time g1 appears in column 1 in pairs, to know its number of neighbors
    ind = gID(ii) == pairs(:,1);      
    if ~isempty(ind)
        gNNeighbors(ii,1) = sum(ind);
    else
        gNNeighbors(ii,1) = 0;
    end
end

gNeighbors = zeros(length(gID), max(gNNeighbors));
for ii = 1:length(gID)
    ind = gID(ii) == pairs(:,1);
    if sum(ind)>0
        gNeighbors(ii, 1:sum(ind)) = reshape(pairs(ind,2),1,[]);
        S.g2{ii} = reshape(pairs(ind,2),1,[]);
    else
        S.g2{ii} = 0;
    end
end


end