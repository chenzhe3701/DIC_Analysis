% chenzhe, 2018-09-07
% Use one_pass_label to generate unique ID map of input matrix. 
% For the uniqueID map, find the ones with grain size smaller than the
% threshold. 
% If they only have one neighbor, then it is a 'hole'. Change its ID to
% that of its only neighbor.

% then return the valid part (i.e., with grain size > threshold) of the input matrix
% min_size, default = 0.25% of input data size

function label = one_pass_fill(raw, min_pct)
%%
if ~exist('min_pct','var')
    min_pct = 0.0025;
end
min_size = sum(raw(:)>0)*min_pct;

% (0) one_pass_label raw map, use negative IDs
label = one_pass_label(raw) + max(raw(:)) + 1;  % make sure label are raw do not have overlap
uniqueLabel = unique(label(:));
uniqueLabel(isnan(uniqueLabel)) = [];
gSize = [];
for ii = 1:length(uniqueLabel)
    gSize(ii) = sum(label(:)==uniqueLabel(ii));
end

% (1) For segments large enough, change the label to the same as in raw.
validLabel = uniqueLabel(gSize>=min_size);
updatedValidLabel = [];
for ii = 1:length(validLabel)
    ind = label == validLabel(ii);
    label(ind) = raw(ind);
    updatedValidLabel =  [updatedValidLabel;unique(raw(ind))];
end
updatedValidLabel = unique(updatedValidLabel);


% (3) find labels corresponding to grains that are too small, set to its
% neighboring ID that are valid

labelTooSmall = uniqueLabel(gSize<min_size);
count = 0;
while sum(ismember(label(:),labelTooSmall))
    count = count + 1;
    for ii = 1:length(labelTooSmall)
        ind = label == labelTooSmall(ii);
        
        ind_nb = imdilate(ind,[0 1 0; 1 1 1; 0 1 0]) - ind;
        nb = label(ind_nb==1);
        nbOK = nb(ismember(nb,updatedValidLabel));
        if ~isempty(nbOK)
            label(ind) = nbOK(1);
            labelTooSmall(ii) = nan;
        end
    end
    labelTooSmall(isnan(labelTooSmall)) = [];
end


end