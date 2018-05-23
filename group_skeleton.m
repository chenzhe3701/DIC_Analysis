% combine consecutive skeletion into a group

function anchor_group = group_skeleton(skl_list)

anchor_group = double(skl_list);
group_number = 0;
pos = 1;
while pos<=length(anchor_group)
   if (skl_list(pos)==1) 
       group_number = group_number + 1;
       while (pos<=length(anchor_group))&&(skl_list(pos)==1)    % sequence is important ...
           anchor_group(pos) = group_number;
           pos = pos + 1;
       end
   end
   pos = pos + 1;
end