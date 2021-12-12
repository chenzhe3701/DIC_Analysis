function initiating_iE_of_each_twin = find_initial_iE_of_twin_in_grain(struCell, ID_current)
% find the iE that each twin initially observed/labeled/detected in this grain.   
% chenzhe, 2019-08-28

% first, find out iE_start
iE_start = 1;
while (isempty(struCell{iE_start}))&&(iE_start<length(struCell))
    iE_start = iE_start + 1;
end

% then, find iE_max
iE_max = iE_start;
while (~isempty(struCell{iE_max}))&&(iE_max<length(struCell))
    iE_max = iE_max + 1;
end
iE_max = iE_max - 1;

% then find iS.  It may happen that ID_current not exist in structure.
iS = find(arrayfun(@(x) x.gID == ID_current, struCell{iE_start}));

% find initiating_iE_of_each_twin
initiating_iE_of_each_twin = inf*ones(1,6);
if ~isempty(iS)
    for iTwin = 1:6
        ie = iE_start;
        while(ie<=iE_max)
            if ~isempty(struCell{ie}(iS).tGb{iTwin})    % ie is the strain level when this twin variant started to be active (having intersecting boundary)
                break;
            else
                ie = ie + 1;
            end
        end
        initiating_iE_of_each_twin(iTwin) = ie;
    end
end

end