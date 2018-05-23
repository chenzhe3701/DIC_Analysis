function [pxl_list_new, ind_list_new, skl_list_new, anchor_label_new]=interp_between_layers(pxl_list_pre, ind_list_pre, skl_list_pre, anchor_label_pre, pxl_list, ind_list, skl_list)

%% for debug
% ind_list_pre = [131, 118, 106, 93, 80, 68, 55, 55, 43, 30, 18, 19, 33, 47, 61, 75, 88, 102, 102, 115, 129, 142, 155, 141, 127, 113, 112, 111, 110, 110, 96, 95, 94, 106, 118];
% skl_list_pre = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1];  
% ind_list = [131, 118, 106, 93, 81, 69, 56, 44, 31, 32, 46, 60, 74, 87, 101, 114, 128, 142, 155, 142, 128, 114, 100, 99, 98, 97, 83, 82, 81, 93, 106, 118];
% pxl_list_pre = rand(size(ind_list_pre));
% pxl_list = rand(size(ind_list));
% skl_list = rand(size(ind_list));
save('for_debug_interp_between_layers.mat','pxl_list_pre', 'ind_list_pre', 'skl_list_pre', 'anchor_label_pre', 'pxl_list', 'ind_list', 'skl_list');

%%
nIntervalsIn = max(group_skeleton(~skl_list));
% anchors_label is the existed skl_label from previous layers
% skl_label is the skl_label for the current layer
skl_label = unique(ind_list(skl_list==1));
% assign group number
anchor_group_pre = group_skeleton(skl_list_pre);
nGroups_pre = max(anchor_group_pre);

% find if it is an Old_Skeleton
skl_list_old = zeros(size(ind_list));
for ii=1:length(skl_list_old)
    if ismember(ind_list(ii), anchor_label_pre)
        skl_list_old(ii) = 1;
    end
end
% assign group number to skl_list Old
anchor_group = group_skeleton(skl_list_old);
nGroups = max(anchor_group);    
% Note that nGroups could < nGroups_pre, which complicates the problem   
% Can nGroups > nGroups_pre? I don't know yet

% Firstly copy, because only need to interp in intervals
pxl_list_new = pxl_list_pre;
ind_list_new = ind_list_pre;

q1_group_num = 1;    % keep track of group number in ind_list, only increase allowed 
for ig = 1:nGroups_pre-1
    p1 = find((skl_list_pre==1)&(anchor_group_pre==ig), 1, 'last');
    p1 = p1 + 1;
    p2 = find((skl_list_pre==1)&(anchor_group_pre==ig+1), 1, 'first');
    p2 = p2 - 1;
    
    % if it is not the final skeleton layer
    if sum(skl_list)~=length(skl_list)     
        % --> possible solution: if q1 and q2 in the same group, then no descendant !!!!! 
        q1 = find((ind_list==ind_list_pre(p1-1))&(anchor_group >= q1_group_num), 1, 'first');
        q1_group_num = anchor_group(q1);
        q1 = find(anchor_group == q1_group_num, 1, 'last');
        q1 = q1 + 1;
        
        q2 = find((ind_list==ind_list_pre(p2+1))&(anchor_group >= q1_group_num), 1, 'last');
        q2_group_num = anchor_group(q2);
        if (q1_group_num==q2_group_num)
            % no descendant
            q1_possible = find((ind_list==ind_list_pre(p1-1))&(anchor_group == q1_group_num)); %, 1, 'first');
            q2_possible = find((ind_list==ind_list_pre(p2+1))&(anchor_group == q1_group_num));

            % find the most similar position
            preMin = inf;
            for ii = 1:length(q1_possible)
                q1_temp = q1_possible(ii);
                for jj = 1:length(q2_possible)
                    q2_temp = q2_possible(jj);
                    if ((q2_temp-q1_temp)>0)&&(q2_temp-q1_temp<=(p2+1)-(p1-1))
                        if min(abs(q1_temp/length(ind_list)-p1/length(ind_list_pre)), abs(q2_temp/length(ind_list)-p2/length(ind_list_pre))) < preMin
                            q1 = q1_temp;
                            q2 = q2_temp;
                        end
                    end
                end
            end
            q1 = q1 + 1;
            q2 = q2 - 1;
            
        else
            q2 = find(anchor_group == q2_group_num, 1, 'first');
            q2 = q2 - 1;
            q1_group_num = q1_group_num + 1;
        end
        
        % perform interp, form [q1 ... q2] to [p1 ... p2]
        if (q2-q1 > p2-p1)
            disp('warning! Interp should not shrink! - 1');
        end
        if (p1==p2)||(q1==q2)
            pos_to_use = q1;
        else
            pos_to_use = interp1(0:q2-q1, q1:q2, linspace(0, q2-q1, p2-p1+1),'nearest');
        end
        ind_list_new(p1:p2) = ind_list(pos_to_use);
        pxl_list_new(p1:p2) = pxl_list(pos_to_use);
    else
        % The Final skeleton layer:
        q1_possible = find(ind_list==ind_list_pre(p1-1));    % an old skeleton point   
        q2_possible = find(ind_list==ind_list_pre(p2+1));    % an old skeleton point        
        % find the most similar position
        preMin = inf;
        for ii = 1:length(q1_possible)
            q1_temp = q1_possible(ii);
            for jj = 1:length(q2_possible)
                q2_temp = q2_possible(jj);
                if ((q2_temp-q1_temp)>0)&&(q2_temp-q1_temp<=(p2+1)-(p1-1))
                   if min(abs(q1_temp/length(ind_list)-p1/length(ind_list_pre)), abs(q2_temp/length(ind_list)-p2/length(ind_list_pre))) < preMin
                      q1 = q1_temp;
                      q2 = q2_temp;
                   end
                end                
            end
        end
        q1 = q1 + 1;
        q2 = q2 - 1;
        
        % perform interp, form [q1 ... q2] to [p1 ... p2]
        if (q2-q1 > p2-p1)
            disp('warning! Interp should not shrink! - 2');
        end
        if (p1==p2)||(q1==q2)
            pos_to_use = q1;
        else
            pos_to_use = interp1(0:q2-q1, q1:q2, linspace(0, q2-q1, p2-p1+1),'nearest');
        end
        ind_list_new(p1:p2) = ind_list(pos_to_use);
        pxl_list_new(p1:p2) = pxl_list(pos_to_use);
    end
    
end

% if has extra data after skeleton
ig = nGroups_pre;
p1 = find(anchor_group_pre==ig,1,'last')+1;
if p1<=length(anchor_group_pre)
    p2 = length(anchor_group_pre);
        
    q1 = find((ind_list==ind_list_pre(p1-1))&(anchor_group >= q1_group_num), 1, 'first');
    currentGroup = anchor_group(q1);
    q1 = find(anchor_group == q1_group_num, 1, 'last');
    q1 = q1 + 1;
    if (q1>length(ind_list))
        disp('Warning! q1 out of range.');
    end
    q2 = find(anchor_group == q1_group_num+1, 1, 'first');
    if ~isempty(q2)
        disp('Warning!, q2 found in another anchor group');
        disp(['p1 group =  ',num2str(anchor_group_pre(p1-1))]);
        disp(['p1, p2, length_list_pre = [',num2str(p1), ',', num2str(p2), ',', num2str(length(anchor_group_pre)),']']);
        disp(['q2 group =  ',num2str(q1_group_num+1)]);
        disp(['q1, q2, length_list = [',num2str(q1), ',', num2str(q2), ',', num2str(length(ind_list)),']']);
        q2 = q2 - 1;
    else
        q2 = length(anchor_group);
    end
    
    % perform interp, form [q1 ... q2] to [p1 ... p2]
    if (q2-q1 > p2-p1)
        disp('warning! Interp should not shrink!');
        disp([p1,p2,q1,q2]);
    end
    if (p1==p2)||(q1==q2)
        pos_to_use = q1;
    else
        try
        pos_to_use = interp1(0:q2-q1, q1:q2, linspace(0, q2-q1, p2-p1+1),'nearest');
        catch
            [p1,p2,q1,q2]
        end
    end
    ind_list_new(p1:p2) = ind_list(pos_to_use);
    pxl_list_new(p1:p2) = pxl_list(pos_to_use);

end

skl_list_new = zeros(size(ind_list_new));
anchor_label_new = union(skl_label, anchor_label_pre);
% update interpolated skeleton list
skl_list_new(ismember(ind_list_new, anchor_label_new)) = 1;

% [debug] check for potential error
nIntervalsOut = max(group_skeleton(~skl_list_new));
if (nIntervalsOut~=nIntervalsIn)
   disp('warning: intervals in and out not consistent : in, out:'); 
   disp([nIntervalsIn,nIntervalsOut]);
end
end


