% function ID_new = hungarian_assign_ID_map(ID, ID_template)
%
% depending on the overlap area, change the id# in [ID] to that in
% [ID_template], and return it as [ID_new]
%
% chenzhe, 2018-05-15
%
% chenzhe, 2018-05-18, fix bug to check both worker and job having match.

function ID_new = hungarian_assign_ID_map(ID, ID_template)

uniqueID = unique(ID(:));
unique_template = unique(ID_template(:));

overlap = zeros(length(uniqueID),length(unique_template));

for ii = 1:length(uniqueID)
   for jj = 1:length(unique_template)
       overlap(ii,jj) = sum( (ID(:)==uniqueID(ii))&(ID_template(:)==unique_template(jj)));
   end
end

[worker, job] = hungarian_assign(max(overlap(:))-overlap);

ID_new = ID;
for ii = 1:length(worker)
    % if this one has a match
    if worker(ii)>0
        id_in_temp = uniqueID(worker(ii));
        % if there is a match, reassign
        if job(ii)~=0
            id_in_template = unique_template(job(ii));
            ID_new(ID_new==id_in_temp) = id_in_template;
        end
    end
end
