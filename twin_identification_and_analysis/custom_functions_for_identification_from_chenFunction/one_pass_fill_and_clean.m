% chenzhe, 2018-09-07
% Use one_pass_label to generate unique ID map of input matrix. 
% For the uniqueID map, find the ones with grain size smaller than the
% threshold. 
% If they only have one neighbor, then it is a 'hole'. Change its ID to
% that of its only neighbor.

% then return the valid part (i.e., with grain size > threshold) of the input matrix
% min_size, default = 0.25% of input data size
%
% chenzhe, 2018-12-17
% add some code to also clean subclusters that are too close to grain boundary by changing cluster
% number to '0'


function label = one_pass_fill_and_clean(raw, size_pct)
%%

% [add new] analyze it again, search for sub-clusters too close to the grain boundary, eliminate  
thisGrain = double(raw~=0); % thisGrain(thisGrain~=1) = -1;    % cannot use 0 to find grain boundary!  Fixed on 2019-01-20.
gbLocal = find_one_boundary_from_ID_matrix(thisGrain);
distMapLocal = bwdist(gbLocal);     % cityDistMapLocal = city_block(gbLocal);     % old code, calculate a cityDistMapLocal 
% try to use a threshold distance
dist_th = quantile(distMapLocal(thisGrain==1), 0.36);     % assume circles 4:5 radius, area ratio is 16:25

if ~exist('size_pct','var')
    size_pct = 0.0025;
end
min_size = sum(raw(:)>0)*size_pct;

% (0) one_pass_label raw map, use negative IDs
label = one_pass_label(raw) + max(raw(:)) + 1;  % make sure label are raw do not have overlap
uniqueLabel = unique(label(:));
uniqueLabel(isnan(uniqueLabel)) = [];
gSize = [];
for ii = 1:length(uniqueLabel)
    gSize(ii) = sum(label(:)==uniqueLabel(ii));
end

% (1) For segments (subclusters) large enough, and distance large enough to the grain boundary, change the label to the same as in raw.
gSizeOK = uniqueLabel(gSize>=min_size);
validLabel = [];
useRaw = false(size(label));
useRaw(thisGrain~=1) = 1; % change cluster number to '0' for areas outside of this grain.  --- Here, outside this code, need to make sure that '0' is used outside of this grain!!!
for ii = 1:length(gSizeOK)
    ind = label == gSizeOK(ii);

    ds = distMapLocal(ind);
    qts = quantile(ds,0.95);
      
    if (qts(end) > dist_th)
        useRaw(ind) = 1;  % useRaw is a logical matrix, 1 indicates the positions where it should be labeled using the cluster number in the 'raw' matrix. 
        validLabel =  [validLabel;unique(label(ind))];  % should record the valid label number in the 'label' matrix. But immediately, this will be changed to the label number in 'raw' matrix, so directly record it in the 'raw' matrix.
    end
end
label(useRaw) = raw(useRaw);
validLabel = unique(validLabel);    % the question is, we need to make sure there is at least one valid? 

if ~any(validLabel==0)
    disp('add 0 to valid label');
    validLabel = [0;validLabel];    % add '0' to validLabel
end

if isempty(validLabel)
   error('empty valid label'); 
end

validLabel_in_updated_mat = unique(label(useRaw(:)==1));
% (3) find labels corresponding to grains that are too small, set to its
% neighboring ID that are valid. (If adjacent to grian boundary, it should
% be set to '0')

nonValidLabel = uniqueLabel(~ismember(uniqueLabel,validLabel));
count = 0;
while sum(ismember(label(:),nonValidLabel))
    count = count + 1;
    for ii = 1:length(nonValidLabel)
        ind = label == nonValidLabel(ii);
        
        ind_nb = imdilate(ind,[0 1 0; 1 1 1; 0 1 0]) - ind;
        nb = label(ind_nb==1);
        nbOK = unique(nb(ismember(nb,validLabel_in_updated_mat)));
        if ~isempty(nbOK)
            label(ind) = nbOK(1);
            nonValidLabel(ii) = nan;
        end
    end
    nonValidLabel(isnan(nonValidLabel)) = [];
end

% label(thisGrain~=1) = 0;

% % [add new] analyze it again, search for sub-clusters too close to the grain boundary, eliminate  
% thisGrain = double(label>0);
% thisGrain(thisGrain~=1) = -1;
% 
% label_cleaned = one_pass_label(label);  % make sure label are raw do not have overlap
% uniqueLabel_cleand = unique(label_cleaned(:));
% uniqueLabel_cleand(isnan(uniqueLabel_cleand)) = [];
% 
% gbLocal = find_one_boundary_from_ID_matrix(thisGrain);
% cityDistMapLocal = city_block(gbLocal);     % calculate a cityDistMapLocal
% gDia = sqrt(sum(thisGrain(:)==1)/pi);
% 
% for ii = 1:length(uniqueLabel_cleand)
%     ind = label_cleaned == uniqueLabel_cleand(ii);
%     ds = cityDistMapLocal(ind);
%     qts = quantile(ds,0.95);
%     % If distance smaller than threshold, eliminate sub-cluster by setting its number to '0'  
%     if (qts(end)<gDia*dist_pct)
%         label_cleaned(ind) = 0;
%     end
% end
% 
% label(label_cleaned==0) = 0;    % if cleaned label has a sub-cluster too close to gb, make it 0

end