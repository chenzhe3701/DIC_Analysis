
% This code was used to fix manual label results.  Adding field 'tGbPtsiE'
%
% chenzhe, 2019-08-15

%% just run this to correct data at iE=2.  Run this and next for data at iE = 2, and 3.  Save by append.
iE = 2
struCell{2}(1).tGbNormal = cell(1,6);
for iS = 1:length(struCell{iE})
    struCell{iE}(iS).tGbPtsiE = cell(1,6);
    for iTwin = 1:6
        % if ~isempty(struCell{iE}(iS).tGb{iTwin})
        for iGb = 1:length(struCell{iE}(iS).tGb{iTwin})     % could be 0, then not run
            for iPt = 1:size(struCell{iE}(iS).tGbPts{iTwin}{iGb},1)
                struCell{iE}(iS).tGbPtsiE{iTwin}{iGb}(iPt,1) = iE;
            end
        end
    end
end

%% 
iE = 3;
struCell{2}(1).tGbNormal = cell(1,6);
struCell{3}(1).tGbNormal = cell(1,6);
for iS = 1:length(struCell{iE})
    struCell{iE}(iS).tGbPtsiE = cell(1,6);
    for iTwin = 1:6
        % if ~isempty(struCell{iE}(iS).tGb{iTwin})
        for iGb = 1:length(struCell{iE}(iS).tGb{iTwin}) 
            for iPt = 1:size(struCell{iE}(iS).tGbPts{iTwin}{iGb},1)
                struCell{iE}(iS).tGbPtsiE{iTwin}{iGb}(iPt,1) = iE;
            end
        end
        
        % if ~isempty(struCell{iE-1}(iS).tGb{iTwin})
        for iGb = 1:length(struCell{iE-1}(iS).tGb{iTwin})
            for iPt = 1:size(struCell{iE-1}(iS).tGbPts{iTwin}{iGb},1)
                struCell{iE}(iS).tGbPtsiE{iTwin}{iGb}(iPt,1) = iE-1;
            end
        end
    end
end
%%
struCell{2}(1).tGbNormal = cell(1,6);
struCell{3}(1).tGbNormal = cell(1,6);

s2=struCell{2}(23); % 23,84
s3=struCell{3}(23);
open s2
open s3