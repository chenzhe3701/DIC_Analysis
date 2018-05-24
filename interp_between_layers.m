% function [pxl_list_new, ind_list_new, skl_list_new, coord_list_new, anchor_label_new]=interp_between_layers(pxl_list_pre, ind_list_pre, skl_list_pre, coord_list_pre, anchor_label_pre, pxl_list, ind_list, skl_list)
function [coord_list_new, anchor_label_new]=interp_between_layers(ind_list_pre, skl_list_pre, coord_list_pre, anchor_label_pre, ind_list, skl_list)

%% for debug
% ind_list_pre = [131, 118, 106, 93, 80, 68, 55, 55, 43, 30, 18, 19, 33, 47, 61, 75, 88, 102, 102, 115, 129, 142, 155, 141, 127, 113, 112, 111, 110, 110, 96, 95, 94, 106, 118];
% skl_list_pre = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1];  
% ind_list = [131, 118, 106, 93, 81, 69, 56, 44, 31, 32, 46, 60, 74, 87, 101, 114, 128, 142, 155, 142, 128, 114, 100, 99, 98, 97, 83, 82, 81, 93, 106, 118];
% pxl_list_pre = rand(size(ind_list_pre));
% pxl_list = rand(size(ind_list));
% skl_list = rand(size(ind_list));
save('for_debug_interp_between_layers.mat', 'ind_list_pre', 'skl_list_pre', 'coord_list_pre', 'anchor_label_pre', 'ind_list', 'skl_list');
% figure; plot(skl_list_pre);
% figure; plot(skl_list);
%%
% [for debugging]
nIntervalsIn = max(group_skeleton(~skl_list));

% skl_label is [a set of (global) index positions] which are skeleton for the current layer, to add to anchor_label 
skl_label = unique(ind_list(skl_list==1));
% anchors_label is the existed skl_label from previous layers, 
% anchor_label_new is updated with info from this layer, output for use for next layer   
anchor_label_new = union(skl_label, anchor_label_pre);

% position_list, for conveniently find/refer/track elements
pos_list = 1:length(ind_list);
pos_list_pre = 1:length(ind_list_pre);

% find mid-point in skl_list and skl_list_pre
mid_ind_list = find_mid_point_in_ind_list(ind_list);
mid_ind_list_pre = find_mid_point_in_ind_list(ind_list_pre);

% first [modify skl_list_pre], so that points with no descendant is also labeled as skeleton.  This is mainly for grouping purpose 
% Make Old_Skeleton list, indicating if current skeleton point already exist in previous layer, which means it is an Old_Skeleton   
skl_list_old = zeros(size(ind_list));
for ii=1:length(skl_list_old)
    if ismember(ind_list(ii), anchor_label_pre)
        skl_list_old(ii) = 1;
    end
end
% Assign group number to skl_list Old
anchor_group = group_skeleton(skl_list_old);
nGroups = max(anchor_group);
% Modify [skl_list_pre] to group segments with no descendants into skeleton  
% if skl_list has [ABC .. mid_point .. CB], then in skl_list_pre: 
% Before mid_piont_pre, anyting A**B, A**C, B**C should be labeled as no descendant.
% After mid_point_pre, anything C****B should be labeled as no descendant.
% To make them no descendant, just change its skl_list_pre to '1'
for ia = 1:mid_ind_list
    if anchor_group(ia)>0
        [~,ib] = find((anchor_group==anchor_group(ia))&(pos_list>ia)&(pos_list<=mid_ind_list),1,'first');
        if ~isempty(ib)
            ic = find((ind_list_pre==ind_list(ia))&(pos_list_pre<=mid_ind_list_pre),1,'first');
            if ~isempty(ic)
                id = find((ind_list_pre==ind_list(ib))&(pos_list_pre>ic)&(pos_list_pre<=mid_ind_list_pre));
                if length(id)>1
                    disp('Warning: should not visit the same point more than twice in half of a traverse');
                    [ia,ib,ic,id]
                end
                if ~isempty(id)
                    skl_list_pre(ic:id) = 1;
                end
            end
        end
    end
end
for ia = mid_ind_list:length(ind_list)
    if anchor_group(ia)>0
        [~,ib] = find((anchor_group==anchor_group(ia))&(pos_list>ia)&(pos_list >= mid_ind_list),1,'first');
        if ~isempty(ib)
            ic = find((ind_list_pre==ind_list(ia))&(pos_list_pre>=mid_ind_list_pre),1,'first');
            if ~isempty(ic)
                id = find((ind_list_pre==ind_list(ib))&(pos_list_pre>ic)&(pos_list_pre>=mid_ind_list_pre));
                if length(id)>1
                    disp('Warning 2: should not visit the same point more than twice in half of a traverse');
                    [ia,ib,ic,id]
                    error('e');
                end
                if ~isempty(id)
                    skl_list_pre(ic:id) = 1;
                end
            end
        end
    end
end


% % modify [skl_list_pre]
% % If skl_list has a pattern [ABCA 000 ...], then in the skl_list_pre, anything between A,B, and A,C should be assigned as skeleton  
% for ia = 1:length(skl_list_old)
%     if anchor_group(ia)>0
%         [~,ib] = find((anchor_group==anchor_group(ia))&(ind_list==ind_list(ia))&(pos_list>pos_list(ia)), 1, 'last');
%         if ~isempty(ib)
%             ic = find(ind_list_pre==ind_list(ia),1,'first');
%             if isempty(ic)
%                 disp('warning: ic not found');
%             end
%             id = find((ind_list_pre==ind_list(ib))&(pos_list_pre>pos_list_pre(ic)), 1, 'last');
%             if ~isempty(id)
%                 skl_list_pre(ic:id) = 1;                
%             end
%         end
%     end
% end
% Also need to do this again reversely, (maybe ...)
% for ia = length(skl_list_old): -1 : 1
%     if anchor_group(ia)>0
%         [~,ib] = find((anchor_group==anchor_group(ia))&(ind_list==ind_list(ia))&(pos_list < pos_list(ia)), 1, 'first');
%         if ~isempty(ib)
%             ic = find(ind_list_pre==ind_list(ia),1,'last');
%             if isempty(ic)
%                 disp('warning: ic not found in reverse search');
%             end
%             id = find((ind_list_pre==ind_list(ib))&(pos_list_pre < pos_list_pre(ic)), 1, 'first');
%             if ~isempty(id)
%                 skl_list_pre(id:ic) = 1;                
%             end
%         end
%     end
% end


% assign group number
anchor_group_pre = group_skeleton(skl_list_pre);
nGroups_pre = max(anchor_group_pre);
% Note that nGroups_pre could be > nGroups, because of no descendant points, 
% which complicates the problem .
% But if skeletons join, nGropus could decrease.
% The above should reduce nGroups_pre (temp_modify) to reduce this problem.
% Can nGroups > nGroups_pre? I don't know yet
if (nGroups < nGroups_pre)
    disp('After re-grouping, nGroups of anchor still reduced!');
end

% % Firstly copy, because only need to interp in intervals
% pxl_list_new = pxl_list_pre;
% ind_list_new = ind_list_pre;
% skl_list_new = skl_list_pre; % this is already the modified skl_list_pre

coord_list_new = zeros(size(ind_list));

% To interp:
% First, find anchor group correspondence
% Then, move to the edge of anchor, and interp the data between anchor groups 
q1_group_num = 1;    % keep track of group number in ind_list, only increase allowed 
q2 = 0;
for ig = 1:nGroups_pre-1
    % [p1 ... p2] is the anchor position in (modified) skl_list_pre
    p1 = find((skl_list_pre==1)&(anchor_group_pre==ig), 1, 'last');
    p2 = find((skl_list_pre==1)&(anchor_group_pre==ig+1), 1, 'first');
    
    if sum(skl_list)~=length(skl_list)
        % If it is not the final skeleton layer
        
        % [q1 ... q2] is the anchor position in skl_list. 
        % find q1:
        % In the same group, but use the last one, (If in different anchor group, use the first one?)  
        q1 = find((ind_list==ind_list_pre(p1))&(anchor_group == q1_group_num), 1, 'last');
        if isempty(q1)
            disp(['Warning: q1 not found in group: ', num2str(ig)]);
        end
        if anchor_group(q1) > q1_group_num
           disp(['Warning: q1 should be found in anchor group: ',num2str(ig),...
               ', but was found in group: ',num2str(anchor_group(q1))]); 
        end
        q1_group_num = anchor_group(q1);
        % find q2:
        q2 = find((ind_list==ind_list_pre(p2))&(anchor_group > q1_group_num), 1, 'first');
        q2 = find((ind_list==ind_list_pre(p2))&(anchor_group >= q1_group_num)&(pos_list>q1), 1, 'first');
        if isempty(q2)
            disp(['Warning: q2 not found in group > than: ', num2str(q1_group_num)]);
        end
        q2_group_num = anchor_group(q2);
        q2 = find((skl_list==1)&(anchor_group==q2_group_num), 1, 'first');
        % error checking, and get interval for interp
%         if (p2-p1) < (q2-q1)
%             disp(['warning, interp should not shrink. [p1,p2]=[',num2str([p1,p2]),'], [q1,q2]=[',num2str([q1,q2]),']']);
%             msg = 'Error due to interp interval shrink';
%             error(msg);
%         end
    else
        % Else, this is the Final skeleton layer
        
        q1_possible = find(ind_list==ind_list_pre(p1));    % an old skeleton point
        q2_possible = find(ind_list==ind_list_pre(p2));    % an old skeleton point
        % find the most similar position
        preMin = inf;
        for ii = 1:length(q1_possible)
            q1_temp = q1_possible(ii);
            for jj = 1:length(q2_possible)
                q2_temp = q2_possible(jj);
                if ((q2_temp-q1_temp)>0)&&(q2_temp-q1_temp<=p2-p1)
                    if min(abs(q1_temp/length(ind_list)-p1/length(ind_list_pre)), abs(q2_temp/length(ind_list)-p2/length(ind_list_pre))) < preMin
                        q1 = q1_temp;
                        q2 = q2_temp;
                    end
                end
            end
        end
        % error checking, and get interval for interp
%         if (p2-p1) < (q2-q1)
%             disp(['warning, interp in Last skeleton layer should not shrink. [p1,p2]=[',num2str([p1,p2]),'], [q1,q2]=[',num2str([q1,q2]),']']);
%         end
    end
    
    % perform the interp
%     p1 = p1 + 1;
%     p2 = p2 - 1;
%     q1 = q1 + 1;
%     q2 = q2 - 1;
    if (p1>p2)
        disp('Error, p1>p2');
    end
%     if (q1 >= q2)
%         disp(['q1=',num2str(q1),' q2=',num2str(q2),' ,interp_pos use q2']);
%         pos_to_use = q2;
%     else
%         pos_to_use = interp1(0:q2-q1, q1:q2, linspace(0, q2-q1, p2-p1+1),'nearest');
%     end
    
    % perform interp, form [q1 ... q2] to [p1 ... p2]
%     ind_list_new(p1:p2) = ind_list(pos_to_use);
%     pxl_list_new(p1:p2) = pxl_list(pos_to_use);
%     skl_list_new(p1:p2) = skl_list(pos_to_use);
    if (q1>q2)
        msg = 'Error, q1 >= q2 when interp position';[q1,q2]
        error(msg);
    elseif (q1==q2)
        disp('when interpolating, q1=q2, so average coord(p1) and coord(p2)');  
        coord_list_new(q1:q2) = (coord_list_pre(p1)+coord_list_pre(p2))/2;
    else
        coord_list_new(q1:q2) = linspace(coord_list_pre(p1), coord_list_pre(p2), q2-q1+1); 
    end
    


    % update q1_group
    q1_group_num = q2_group_num;

end

% if has extra data after skeleton
ig = nGroups_pre;
p1 = find((skl_list_pre==1)&(anchor_group_pre==ig),1,'last');

if p1 < length(anchor_group_pre)
%     disp('Extra data after last anchor_group');
    p2 = length(anchor_group_pre);
    
    % based on previous q2, find q1 position (should change to last one in the same group) 
    q1 = find((ind_list==ind_list_pre(p1))&(pos_list>q2), 1, 'first');
    if isempty(q1)
        disp('Warning! q1 not found .');
    end
    q1_group_num = anchor_group(q1);
    q1 = find((ind_list==ind_list_pre(p1))&(anchor_group == q1_group_num), 1, 'last');
    
    q2 = length(anchor_group);

    % perform interp, form [q1 ... q2] to [p1 ... p2]
%     p1 = p1 + 1;
%     q1 = q1 + 1;
    if (p1>p2)
        disp('Error, p1>p2');
    end
%     if (q2-q1 > p2-p1)
%         disp('warning! Interp after skeleton should not shrink!');
%         disp([p1,p2,q1,q2]);
%     end
    if (q1 >= q2)
        disp(['q1=',num2str(q1),' q2=',num2str(q2),' ,interp_pos use q2']);
        pos_to_use = q2;
    else
        pos_to_use = interp1(0:q2-q1, q1:q2, linspace(0, q2-q1, p2-p1+1),'nearest');
    end
    
    
    if (q1>q2)
        msg = 'Check error, here q1 should not larger than q2';
        error(msg);
    elseif (q1==q2)
        disp('when interpolating data after all skeleton, q1=q2, so average coord(p1) and coord(p2)');
        disp(['In this condition, it should be p1=p2, here [p1,p2]=[',num2str([p1,p2]),']']);
        coord_list_new(q1:q2) = (coord_list_pre(p1)+coord_list_pre(p2))/2;
    else
        coord_list_new(q1:q2) = linspace(coord_list_pre(p1), coord_list_pre(p2), q2-q1+1);
    end
    
        
        
    % perform interp, form [q1 ... q2] to [p1 ... p2]
%     ind_list_new(p1:p2) = ind_list(pos_to_use);
%     pxl_list_new(p1:p2) = pxl_list(pos_to_use);
%     skl_list_new(p1:p2) = skl_list(pos_to_use);

end

% Now need to treat the unassigned 0's in coord_list_new
group = group_skeleton(~logical(coord_list_new));
if coord_list_new(1)==0
    coord_list_new(1) = coord_list_pre(1);
end
if (coord_list_new(end)~=0)&&(coord_list_new(end)~=coord_list_pre(end))
    disp(['warning, coord_list_new(end) = ',num2str(coord_list_new(end))]);
    disp(['which != coord_list_pre(end) = ',num2str(coord_list_pre(end))]);
end
if coord_list_new(end)==0
    coord_list_new(end) = coord_list_pre(end);
end
for ig = 1:max(group(:))
   p1 = find(group==ig,1,'first');
   p2 = find(group==ig,1,'last');
   if(p1>1)
       p1 = p1-1;
   end
   if(p2<length(group))
       p2 = p2 + 1;
   end
   
   coord_list_new(p1:p2) = linspace(coord_list_new(p1), coord_list_new(p2), p2-p1+1);
end

% skl_list_new2 = zeros(size(ind_list_new));
% % Update interpolated skeleton list
% skl_list_new2(ismember(ind_list_new, anchor_label_new)) = 1;
% % double chek
% if sum(skl_list_new2-skl_list_new)
%     disp('Warning: skl_list_new by two methods do not match');
% end
% % [debug] check for potential error
% nIntervalsOut = max(group_skeleton(~skl_list_new));
% if (nIntervalsOut~=nIntervalsIn)
%    disp('warning if chosen to replace skl_list with skl_list_new, because intervals in and out not consistent : in, out:'); 
%    disp([nIntervalsIn,nIntervalsOut]);
%    msg = 'Error due to in-consistent in and out intervals';
%    error(msg);
% end

end


