function mid_point = find_mid_point_in_ind_list(ind_list)

pos_list = 1:length(ind_list);

p = length(ind_list);
mid_point = p;
for p1 = 1:length(ind_list)
    p = find(ind_list==ind_list(p1));
    if length(p)>=2
        break;
    end    
end
if length(p)>1
    p1 = p(1);
    p2 = p(2);
    while(p1<p2)
       p1 = p1 + 1;
       p2_trial = find((ind_list==ind_list(p1))&(pos_list>p1),1,'last');
       if ~isempty(p2_trial)
           p2 = p2_trial;
           mid_point = p2;
       end
    end
end
if ind_list(mid_point-1)~=ind_list(mid_point)
    mid_point = mid_point - 1;
end


end

% safer version
% p = length(skl_list);
% mid_point = p;
% for p1 = 1:length(skl_list)
%    if (skl_list(p1)==1)
%        p = find(ind_list==ind_list(p1));
%        if length(p)>=2
%            break;
%        end
%    end       
% end
% if length(p)>1
%     p1 = p(1);
%     p2 = p(2);
%     while(p1<p2)
%        p1 = p1 + 1;
%        if (skl_list(p1)==1)
%            p2_trial = find((ind_list==ind_list(p1))&(pos_list>p1),1,'last');
%            if ~isempty(p2_trial)
%                p2 = p2_trial;
%                mid_point = p2;
%            end
%        end
%     end
% end