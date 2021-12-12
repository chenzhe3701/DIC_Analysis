% input the iE and iCluster number, and struCell
% output the tracked iE_list and iC_list for the cluster
%
% chenzhe, 2019-09-17, specifically used for cluster tracking codes.

function [iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE, iCluster)

iE_list = iE;
iC_list = iCluster;


% search to earlier strain (or try not to...)
while 0 ~= struCell{iE_list(1)}(iS).preCluster(iC_list(1))
    iC_list = [struCell{iE_list(1)}(iS).preCluster(iC_list(1)), iC_list];
    iE_list = [iE_list(1)-1, iE_list];
end
% search to later strain
while 0 ~= struCell{iE_list(end)}(iS).postCluster(iC_list(end))
    iC_list = [iC_list,struCell{iE_list(end)}(iS).postCluster(iC_list(end))];
    iE_list = [iE_list, iE_list(end)+1];
end

end