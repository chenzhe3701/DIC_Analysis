% function ID_new = hungarian_assign_ID_map(ID, ID_template)
%
% depending on the overlap area, change the id# in [ID] to that in
% [ID_template], and return it as [ID_new]
%
% chenzhe, 2018-05-15
%
% chenzhe, 2018-05-18, fix bug to check both worker and job having match. 
% chenzhe, 2018-06-14, modify to make it faster, using accumarray
% chenzhe, 2018-06-27, add output, ID_new_unbalance can try to assign as
% many as possible. (ID_new is one-to-one match).
% So, if clean EBSD, ID_new_unbalance can be used to assign euler angle.
% And ID_new can be used to assign ID.

function [ID_new, ID_new_unbalance] = hungarian_assign_ID_map(ID, ID_template)

% should ignore ID = 0. First if there are nans, convert to 0
ID(isnan(ID)) = 0;
ID_template(isnan(ID_template)) = 0;

[uniqueID,~,m1]=unique(ID);
[unique_template,~,m2]=unique(ID_template);

overlap = accumarray([m1,m2],1);

% % This is easy to understand, but very slow.
% % It is even slower if do not find all the pairs first.
% 
% uniqueID = unique(ID(:));
% unique_template = unique(ID_template(:));
% 
% pair = unique([ID(:),ID_template(:)],'rows');
% pair(isnan(sum(pair,2)),:) = [];
% 
% overlap = zeros(length(uniqueID),length(unique_template));
% 
% for ii = 1:length(uniqueID)
%    for jj = 1:length(unique_template)
%        if ismember([uniqueID(ii),unique_template(jj)],pair,'rows')
%         overlap(ii,jj) = sum( (ID(:)==uniqueID(ii))&(ID_template(:)==unique_template(jj)));
%        end
%    end
% end

[worker, job, worker_full, job_full] = hungarian_assign(max(overlap(:))-overlap);
save('worker_job.mat','worker', 'job', 'worker_full', 'job_full');

% worker = 1:size(overlap,1);
% job = munkres(max(overlap(:))-overlap);

% (1) one-to-one match
ID_new = ID;
for ii = 1:length(worker)
    % if this one has a match
    if worker(ii)>0
        id_in_temp = uniqueID(worker(ii));
        % if there is a match, reassign
        if job(ii)~=0
            id_in_template = unique_template(job(ii));
            ID_new(ID==id_in_temp) = id_in_template;
        end
    end
end

% (2) all matched
ID_new_unbalance = ID;
for ii = 1:length(worker_full)
    % if this one has a match
    if worker_full(ii)>0
        id_in_temp = uniqueID(worker_full(ii));
        % if there is a match, reassign
        if job_full(ii)~=0
            id_in_template = unique_template(job_full(ii));
            ID_new_unbalance(ID==id_in_temp) = id_in_template;
        end
    end
end

    
disp('');
