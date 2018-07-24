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

function [ID_new, id_link_additional] = hungarian_assign_ID_map(ID_temp, ID_target)

% should ignore ID = 0. First if there are nans, convert to 0
ID_temp(isnan(ID_temp)) = 0;
ID_target(isnan(ID_target)) = 0;

[uniqueID_temp,~,m1]=unique(ID_temp);
[uniqueID_target,~,m2]=unique(ID_target);

overlap = accumarray([m1,m2],1);

% % This is easy to understand, but very slow.
% % It is even slower if do not find all the pairs first.
% 
% uniqueID_temp = unique(ID_temp(:));
% uniqueID_target = unique(ID_target(:));
% 
% pair = unique([ID_temp(:),ID_target(:)],'rows');
% pair(isnan(sum(pair,2)),:) = [];
% 
% overlap = zeros(length(uniqueID_temp),length(uniqueID_target));
% 
% for ii = 1:length(uniqueID_temp)
%    for jj = 1:length(uniqueID_target)
%        if ismember([uniqueID_temp(ii),uniqueID_target(jj)],pair,'rows')
%         overlap(ii,jj) = sum( (ID_temp(:)==uniqueID_temp(ii))&(ID_target(:)==uniqueID_target(jj)));
%        end
%    end
% end

[worker, job, worker_full, job_full] = hungarian_assign(max(overlap(:))-overlap);
save('worker_job.mat','worker', 'job', 'worker_full', 'job_full');

% Find out the unbalanced, additional match. i.e., after one-to-one match,
% assign as many matches as possible
ind_additional_match = ~ismember(worker_full,worker);
worker_a = worker_full(ind_additional_match);
job_a = job_full(ind_additional_match);


% worker = 1:size(overlap,1);
% job = munkres(max(overlap(:))-overlap);

% (1) one-to-one match
ID_new = ID_temp;
for ii = 1:length(worker)
    % if this one has a match
    if worker(ii)>0
        id_in_temp = uniqueID_temp(worker(ii));
        % if there is a match, reassign
        if job(ii)~=0
            id_in_target = uniqueID_target(job(ii));
            ID_new(ID_temp==id_in_temp) = id_in_target;
        end
    end
end


% (2) after one-to-one match, assign additional matches
id_additional = max(uniqueID_target) + 1;
id_link_additional = [];

for ii = 1:length(worker_a)
    % if this one has a match
    if worker_a(ii)>0
        id_in_temp = uniqueID_temp(worker_a(ii));
        % if there is a match, reassign
        if job_a(ii)~=0
            id_in_target = uniqueID_target(job_a(ii));
            % Make additional assignment, and record the linkage for additional assignment
            ID_new(ID_temp==id_in_temp) = id_additional;
            id_link_additional = [id_link_additional; id_additional, id_in_target];
            % after using, increment id_additional for the next assignment
            id_additional = id_additional + 1;
            
            % method-2, directly find the max in overlap matrix. This just confirms that it is the result we want. 
            % [val,ind] = max(overlap(worker_a(ii),:));
            % disp([ind,job_a(ii),ind-job_a(ii),val]);
        end
    end
end


% % (3) all matched. This is old code.
% ID_new_unbalance = ID_temp;
% for ii = 1:length(worker_full)
%     % if this one has a match
%     if worker_full(ii)>0
%         id_in_temp = uniqueID_temp(worker_full(ii));
%         % if there is a match, reassign
%         if job_full(ii)~=0
%             id_in_target = uniqueID_target(job_full(ii));
%             ID_new_unbalance(ID_temp==id_in_temp) = id_in_target;
%         end
%     end
% end

    
disp('');
