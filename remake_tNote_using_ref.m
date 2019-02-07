
% here I need a code to quickly compare the previously corrected result to
% the new result.
% (1) If it was previously (i) identified or (ii) corrected --> check
% stru.cTrueTwin(iC x iTwin)
% but not identified now, try to make the corresponding tNote as if this is
% done interactively.
% chenzhe, 2018-12-30.

predata = '20190204_0607_relabeled_result.mat';
predata = matfile(predata);

%% 
struCell_ref = predata.struCell;
% trueTwinMapCell_ref = predata.trueTwinMapCell;
tNote = zeros(1,9);
for iE = iE_start:iE_stop
    for iS = 1:length(struCell{iE})
        for iC = 1:size(struCell{iE}(iS).cTrueTwin,1)
            if any(struCell_ref{iE}(iS).cTrueTwin(iC,:) - struCell{iE}(iS).cTrueTwin(iC,:))
                v_input = [struCell{iE}(iS).gID, iE, iC, struCell_ref{iE}(iS).cTrueTwin(iC,:)];
                tNote = [tNote;v_input];
            end
        end
    end
end
tNote = unique(tNote,'rows');
size(tNote)