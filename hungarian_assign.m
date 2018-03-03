% hungarian_assign algorithm, naive
% m = cost_of(worker,job)
% output is the job that each worker should do.  
% Workers were sorted ascending, so worker is omitted in default output.
%
% chenzhe, 2018-03-02

function [worker, job, worker_full, job_full] = hungarian_assign(M)

if ~exist('M','var')
    n = 5;
    disp('worst case example for 10x10 matrix: ');
    M = [1:n]'.*[1:n]
end

M0 = M;
[nR,nC] = size(M);
if nR<nC
    M(nR+1:nC,:) = 0;
elseif nR>nC
    M(:,nC+1:nR) = 0;
end

tf = false;
nLoops = 0;
while ~tf
    nLoops = nLoops + 1;
    % (1) try to assign by optimizing each row.
    M = M - min(M,[],2);
    [tf, worker, job] = check_assignment(M);
    
    % (2) If does not work, try to assign by optimizing each col.
    if tf
        break;
    else
        M = M - min(M,[],1);
        [tf, worker, job] = check_assignment(M);
    end
    
    % (3) if not solved, keep going.
    % Cover lines with 0's with as few rows/cols as possible. One way of doing it:
    if tf
        break;
    else
        assigned = false(size(M));
        crossed = false(size(M));
        marked = false(size(M));
        lined = false(size(M));
        % [1] First, Assign as many tasks as possible. (Already did in check_assignment) 
        % and cross all other 0s within the same row and col
        for ii = 1:length(worker)
            iR = worker(ii);
            iC = job(ii);
            assigned(iR,iC) = true;
            crossed(iR,(0==M(iR,:))&(~assigned(iR,:)))=true;
            crossed((0==M(:,iC))&(~assigned(:,iC)),iC)=true;
        end
        
        % [2] Now, to the drawing part
        marked((0==sum(assigned,2)),:) = true;  % mark all rows having no assigment
        
        cols = sum((0==M)&(marked),1)>0;    % mark all (unmarked) columns having zeros in 'newly marked' rows (so, used 'marked' in index)
        cols_old = [];
        while ~isequal(cols,cols_old)
            marked(:,cols) = true; 
            marked(sum((assigned)&(marked),2)>0,:) = true; % mark all rows having assignments in newly marked columns
            cols_old = cols;
            cols = sum((0==M)&(marked),1)>0;
        end
        % Repeat for all new rows with 0s, until no more cols need to be marked 
        
        % draw lines through all makred columns and unmarked rows
        lined(:,size(marked,1)==sum(marked,1)) = true;
        lined(size(marked,2)~=sum(marked,2),:) = true;
        
        % From the elements that are left, find the lowest value.
        % Subtract this from every unmarked element and add it to every element covered by two lines.
        lowest = min(M(~lined));
        M(~lined) = M(~lined) - lowest;
        M(size(lined,2)==sum(lined,2), size(lined,1)==sum(lined,1)) = ...
            M(size(lined,2)==sum(lined,2), size(lined,1)==sum(lined,1)) + lowest;
    end
    
end

worker_full = worker;
job_full = job;

% if have virtual worker, the job can be done by some real worker
% if nR<nC
for ii = 1:length(worker_full)
    if worker_full(ii)>nR
        worker(ii) = 0;
        [~,worker_full(ii)] = min(M0(1:nR,job_full(ii)));
    end
end

% if have virtual job, the worker can do some real job
% (if nR>nC)
for ii = 1:length(worker_full)
    if(job_full(ii)>nC)
        job(ii) = 0;
        [~,job_full(ii)] = min(M0(worker_full(ii),1:nC));
    end
end


end


function [tf, worker, job] = check_assignment(M)

if ~isempty(M)
    % try to reduce one row and one col
    zeros_in_col = sum(0==M,1);
    if any(zeros_in_col)
        zeros_in_col(0==zeros_in_col) = inf;    % prevent from finding
        [~,indc] = min(zeros_in_col);     % find the col with least 0's
    else
        worker = [];
        job = [];
        tf = false;
        return;
    end
    indrs = find(M(:,indc)==0);    % find a row with least 0's in the col
    [~,ind] = min(sum(0==M(indrs,:),2));
    indr = indrs(ind);  % !!!
    M(:,indc) = [];  % mark col
    M(indr,:) = [];  % mark row
    
    if isempty(indr)
        tf = false;
        return;
    end
    [tf, worker, job] = check_assignment(M);  % find assignment in sub matrix
    % for index larger than the marked row/column, add index by 1
    worker(worker>=indr) = worker(worker>=indr) + 1;
    job(job>=indc) = job(job>=indc) + 1;
    % append the assignment in sub matrix to the current assignment in this iteration
    worker = [indr,worker];
    job = [indc,job];
    % sort ascend
    wj = sortrows([worker(:),job(:)]);
    wj = wj';
    worker = wj(1,:);
    job = wj(2,:);
else
    worker = [];
    job = [];
    tf = true;
end

end

