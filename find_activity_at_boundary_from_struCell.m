function [iE_iTwin_dToTriple_list, has_grain_data] = find_activity_at_boundary_from_struCell(struCell, iE, ID_current, gb, xy_triple)
% function to return activities in this grain (ID_current) at this boundary
% (gb), up to this iE (iE), and dist of intersection to triple point: 
% [iEs, iTwins, dist_to_triple, triple_point_x, triple_point_y].  
% If no activity, return inf.
%
% chenzhe, 2019-08-27 modify

iS = find(arrayfun(@(x) x.gID == ID_current, struCell{iE}));
iE_iTwin_dToTriple_list = [];
if isempty(iS)
    has_grain_data = false;
%     disp(['Grain ',num2str(ID_current), ' data not found']);
else
    has_grain_data = true;
    for iTwin = 1:6
        if ismember(gb, struCell{iE}(iS).tGb{iTwin})
            iGb = find(gb == struCell{iE}(iS).tGb{iTwin});  % find if iTwin intersect the target gb
            tGbPts = struCell{iE}(iS).tGbPts{iTwin}{iGb};      % find the pts
            tGbPtsiE = struCell{iE}(iS).tGbPtsiE{iTwin}{iGb};  % find the iEs
            
            % for each intersecting point (each correspond to an iE)
            for ii = 1:length(tGbPtsiE)
                [dToTriple, ind] = min(pdist2(tGbPts(ii,:), xy_triple), [], 2);
                iE_iTwin_dToTriple_list = [iE_iTwin_dToTriple_list; tGbPtsiE(ii), iTwin, dToTriple, xy_triple(ind,:)];
            end
        end
    end
    iE_iTwin_dToTriple_list = sortrows(iE_iTwin_dToTriple_list);
end

if isempty(iE_iTwin_dToTriple_list)
    iE_iTwin_dToTriple_list = [inf, inf, inf, inf, inf];
end

end