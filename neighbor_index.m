% ind_nb = nb_idx(nR,nC,ind,ind_opt)
% Given the index = ind, or sub=[ind,ind_opt] of a point
% return its neighbor's index
%
% chenzhe, 2018-05-20

function ind_nb = neighbor_index(nR,nC,ind,ind_opt)

if exist('ind_opt','var')
    iR = ind;
    iC = ind_opt;
else
    [iR,iC] = ind2sub([nR,nC],ind);
end

% A point (ind) at [iR,iC] has <= 8 neighbors, the ind differences are:
% (-nR-1)*(iC>1)*(iR>1)     (-1)*(iR>1)     (+nR-1)*(iC<nC)*(iR>1)
% (-nR)*(iC>1)              0               (+nR)*(iC<nC)
% (-nR+1)*(iC>1)*(iR<nR)    (+1)*(iR<nR)    (+nR+1)*(iC<nC)*(iR<nR)

% difference in index, counter clockwise
ind_diff = [(-nR-1)*(iC>1)*(iR>1), (-nR)*(iC>1), (-nR+1)*(iC>1)*(iR<nR), (+1)*(iR<nR), ...
    (+nR+1)*(iC<nC)*(iR<nR), (+nR)*(iC<nC), (+nR-1)*(iC<nC)*(iR>1), (-1)*(iR>1)];
% eliminate zeros, and redundancy
ind_diff = ind_diff(ind_diff~=0);

ind_nb = ind + ind_diff;


end
