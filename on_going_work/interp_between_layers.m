% function [pxl_list_new, ind_list_new, skl_list_new, coord_list_new, anchor_label_new]=interp_between_layers(pxl_list_pre, ind_list_pre, skl_list_pre, coord_list_pre, anchor_label_pre, pxl_list, ind_list, skl_list)
function [coord_list_new]=interp_between_layers(ind_list_pre, skl_list_pre, coord_list_pre, ind_list, skl_list)

%% for debug
% ind_list_pre = [131, 118, 106, 93, 80, 68, 55, 55, 43, 30, 18, 19, 33, 47, 61, 75, 88, 102, 102, 115, 129, 142, 155, 141, 127, 113, 112, 111, 110, 110, 96, 95, 94, 106, 118];
% skl_list_pre = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1];  
% ind_list = [131, 118, 106, 93, 81, 69, 56, 44, 31, 32, 46, 60, 74, 87, 101, 114, 128, 142, 155, 142, 128, 114, 100, 99, 98, 97, 83, 82, 81, 93, 106, 118];
% pxl_list_pre = rand(size(ind_list_pre));
% pxl_list = rand(size(ind_list));
% skl_list = rand(size(ind_list));
save('for_debug_interp_between_layers.mat', 'ind_list_pre', 'skl_list_pre', 'coord_list_pre', 'ind_list','skl_list');
% figure; plot(skl_list_pre);
% figure; plot(skl_list);

%%
% position_list, for conveniently find/refer/track elements
pos_list = 1:length(ind_list);
pos_list_pre = 1:length(ind_list_pre);
coord_list_new = zeros(size(ind_list));

p1 = 0; % left
p2 = 1; % right
q1 = 0;
q2 = 1;

while(p1<length(ind_list_pre))
    p1 = find((skl_list_pre==1)&(pos_list_pre >= p2),1,'first');
    if ~isempty(p1)
        q1 = find((ind_list==ind_list_pre(p1))&(pos_list >= q2),1,'first');
        if isempty(q1)
            error('Error 1, found p1 but cannot find q1');
        end
        
        p2 = find((skl_list_pre==1)&(pos_list_pre > p1), 1, 'first');
        if ~isempty(p2)
           q2 = find((ind_list==ind_list_pre(p2))&(pos_list > q1),1,'first');
            if isempty(q2)
                error('Error 2, found p2 but cannot find q2');
            
            end
            % special treatment just try
            if ((q2-q1)+10<(p2-p1))&&(length(find((ind_list==ind_list_pre(p2))&(pos_list > q1)))>1)
               q2_temp =  find((ind_list==ind_list_pre(p2))&(pos_list > q1));
               [~,ind] = min(abs(q2_temp/length(ind_list)-p2/length(ind_list_pre)));
               q2 = q2_temp(ind);
            end
                
        end
    end
    
    if isempty(p2)&&(p1<length(ind_list_pre))
        p2 = length(ind_list_pre);
        q2 = length(ind_list);
        if ((p1==p2)||(q1==q2))
            error('condition not considered');
        end
    end
   
    % interp
    coord_list_new(q1:q2) = linspace(coord_list_pre(p1), coord_list_pre(p2), q2-q1+1);
    p1 = p2;
    q1 = q2;
end


end


