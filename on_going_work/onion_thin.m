% [pxl_M, ind_M, skl_M] = onion_thin(A, ind_start_input)
% A = input image
% ind_start_point is the index to start traversal.
% chenzhe, 2018-05-19



function [pxl_M, ind_M, skl_M] = onion_thin(A, ind_start_input)

save('for_debug_onion_thin.mat','A');
% For debug.
% clear; clc; close all;

% Use grain ID for debugging:
realGrain = 0;
if realGrain
    ID = load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt\WE43_T6_C1_s_all_grain_1144_local_map.mat');
    ID = ID.data;
    ID = ID(5).ID_local;
    ID = ID==1144;
    A = ID;
    ind_start_input = find(A(:)==1,1,'first');
end

if ~exist('A','var')
    load('A.mat','A');
    ind_start_input = 79;
end
if ~exist('ind_start_input','var')
    ind_start_input = find(A(:)>0,1,'first');
end

[nR,nC] = size(A);
ind_mat = reshape(1:nR*nC, nR, nC);
% find neighbors index
nb_ind_cell = arrayfun(@(x) {neighbor_index(nR,nC,x)}, ind_mat);
% The final skeleton
skeleton = thin(A,inf);

% [debug] make a copy of the original A, make pixel value with index number
% A_raw = double(A).*ind_mat;  

A_raw = double(A);  
A = logical(A);

pxl_cell = [];
ind_cell = [];
skl_cell = [];

% start at
ind = ind_start_input;
reAligned = false;

%% Enter level-1 loop, for each layer
all_skeleton = 0;
iLoop = 1;
while (~all_skeleton)&&(iLoop<=inf)
    
    pxl_list = [];
    ind_list = [];
    skl_list = [];
    coord_list = [];
    
    % Thin by one layer
    A_thin_once = thin(A);
    del_this_loop = logical(A) - A_thin_once;
    
    % initialize setting for this layer
    ind_layer_start = ind;
    ind_layer_2nd_pt = 0;
    
    % If not reAligned, need to assume an old position
    if ~reAligned
        nb_ind_list = nb_ind_cell{ind_layer_start};
        ind_previous = nb_ind_list(end);
    end
    
    %% Enter level-2 loop, traverse each pixel on the boundary of this layer
    while true
        % Find next neighbor:
        % get neighbor index list
        nb_ind_list = nb_ind_cell{ind};
        % circshift so that the last step's ind moves to the end of list
        shift_amount = find(nb_ind_list == ind_previous);
        if isempty(shift_amount)
            shift_amount = 0;
        end
        nb_ind_list = circshift(nb_ind_list, -shift_amount, 2);
        % Find the first foreground pixel, which is the next pixel on GB.
        ind_next = nb_ind_list(find(A(nb_ind_list)==1, 1, 'first'));
        % if it is reduced to a single point
        if isempty(ind_next)
            ind_next = ind;
        end
        
        % stop traversal this layer
        if (ind==ind_layer_start)&&(ind_next==ind_layer_2nd_pt)
            break;  % stop at [ind_previous, ind, ind_next];
        end
        
        % Assign the next step of the layer start position
        if (ind_layer_2nd_pt==0)
            ind_layer_2nd_pt = ind_next;
        end
        
        % record necessary value
        pxl_list = [pxl_list, A_raw(ind)];
        ind_list = [ind_list, ind]; 
        skl_list = [skl_list, skeleton(ind)];
        coord_list = [coord_list, length(pxl_list)-1];
        
        % update ind to next position
        ind_previous = ind;
        ind = ind_next;
        
        % if (realigned==False)&&(find(skl_list)), realign, then realigned=true
    end
    
    % Update 'A' after this layer was traversed
    % Break if only skeleton is left
    if sum(A(:)-A_thin_once(:))
        A = A_thin_once;
    else
        disp(['Finished traversing after Loop: ', num2str(iLoop)]);
        all_skeleton = 1;
    end
    
    %     % [help debug, show the path]
    %     path = zeros(size(A));
    %     path(ind_list) = 1:length(ind_list);
    %     figure;imagesc(path);title('path');
    %     figure;imagesc(del_this_loop);title('delete this loop');
    %     figure;imagesc(A);title('newA');
    
    % if not realigned, do these:
    if ~reAligned
        % interpolate pxl_val_list, ind_list, skl_list
        if iLoop==1
           list_length = length(ind_list);  % 1st layer does not need to interp, but it can define the width of the new map (list_length) 
           coord_range = list_length - 1;
           coord_list = 0:coord_range;
        else
            % if only have one point in this layer, make it 2 points to interp 
            if length(ind_list)==1
                pxl_list = [pxl_list,pxl_list];
                ind_list = [ind_list,ind_list];
                skl_list = [skl_list,skl_list];
                coord_list = [0,1];
            end
            
            % scale the coordinate position to the range [0, coord_range]
            coord_list = coord_list/range(coord_list)*coord_range;
            
            % pxl_list = interp1(1:length(pxl_list), pxl_list, linspace(1,length(pxl_list),list_length),'nearest');
            % ind_list = interp1(1:length(ind_list), ind_list, linspace(1,length(ind_list),list_length),'nearest');
            % skl_list = interp1(1:length(skl_list), skl_list, linspace(1,length(skl_list),list_length),'nearest');
        end
        
        % find starting point of next layer
        ind = nb_ind_list(find(A(nb_ind_list)==1, 1, 'first'));
        
        % try to see if can re-align
        align_shift = find(skl_list, 1, 'first')-1;
        
        % It can be re-aligned
        if ~isempty(align_shift)
            reAligned = true;
            pxl_list = circshift(pxl_list, -align_shift, 2);
            ind_list = circshift(ind_list, -align_shift, 2);
            skl_list = circshift(skl_list, -align_shift, 2);
            % coord_list stay the same 
            
            layer_at_realign = iLoop;
            layer_range_at_realign = length(ind_list)-1;
            ind = ind_list(1);
            ind_previous = ind_list(end);
            anchor_label_cell = [];
            anchor_label_cell{iLoop} = unique(ind_list(skl_list==1));
        end
        % copy the list into cell. If realigned, then copy the realigned version, and record info: [align_shift, layer_at_realign]     
        pxl_cell{iLoop} = pxl_list;
        ind_cell{iLoop} = ind_list;
        skl_cell{iLoop} = skl_list;
        coord_cell{iLoop} = coord_list;
    else
        % align by considering anchor points
        % Use last loops aligned ind_list as template, find all the anchor positions in index in this list, make sure the last one is considered as anchor  
        % Only interp between different GROUPs of anchor points.
        % --> If two anchor points are in the same group in step-i, they will remain in the same group in step i+1.
        % It is possible that some anchor points in a group in step i can show twice in step i+1.
        % But if we interp using the largest range, because anchors deccendants will be eliminated, so it is still OK to just interp.  
        % I think by finding the 'last' of p2, it can be achieved that two anchors are from different 'groups'  
        % example-1: nIntervals decrease, some doesn't have descendants any more   
        % layer-1: [14, 15, 16, ... ...]
        % layer-2: [14, 16, 14, ...]
        % example-2: ?
        % layer-1: [226, 225, ..., ..., ..., ..., ..., 715,225]
        % layer-2: [226, 225, 715, ..., 265, 715, 225]
        % example - 3:
        % layer-1: [131,118,106,93, ..., ..., ..., ..., ..., 142, 155, ..., (x10), ..., 106,118]  
        % layer-2: [131,118,106,93, ..., 142, 155, 142, ..., (x10), ..., 106, 118]
        ind_list_old = ind_list;
        [coord_list] ...
            = interp_between_layers(ind_cell{iLoop-1}, skl_cell{iLoop-1}, coord_cell{iLoop-1},ind_list,skl_list);

        % debug, check when some points are missing
        if length(unique(ind_list))<length(unique(ind_list_old))
            close;
            disp(['found point missing at loop: ',num2str(iLoop),' ----------']);
            save('debugdata_point_missing','iLoop','ind_cell','skl_cell','ind_list','ind_list_old');
            myplot(A);
        end
            
        % copy the list into cell
        pxl_cell{iLoop} = pxl_list;
        ind_cell{iLoop} = ind_list;
        skl_cell{iLoop} = skl_list;
        coord_cell{iLoop} = coord_list;
    end
    
    % increment loop
    iLoop = iLoop + 1;
end

%% check range
pos = [];
for ii = 1:length(coord_cell)
    pos = [pos; coord_cell{ii}(1), coord_cell{ii}(end)];     
end
if (pos(1,1)~=max(pos(:,1)))||(pos(1,2)~=max(pos(:,2)))
    disp('warning! finished thinning but inner layer seems larger !');
end

%% align
for ii = 1:length(coord_cell)
    u_target = 0:coord_range;
    
    F = griddedInterpolant(coord_cell{ii},pxl_cell{ii},'nearest');
    pxl_cell{ii} = F(u_target);
    F = griddedInterpolant(coord_cell{ii},ind_cell{ii},'nearest');
    ind_cell{ii} = F(u_target);
    F = griddedInterpolant(coord_cell{ii},skl_cell{ii},'nearest');
    skl_cell{ii} = F(u_target);
end



%% make into matrix
make_matrix = 1;    % for debug purpose. Otherwise just return the cell.
if make_matrix
    pxl_M = cell2mat(reshape(pxl_cell,[],1));
    ind_M = cell2mat(reshape(ind_cell,[],1));
    skl_M = cell2mat(reshape(skl_cell,[],1));
    
    % shift back
    align_shift = round(align_shift/layer_range_at_realign*coord_range);
    pxl_M(layer_at_realign:end,:) = circshift(pxl_M(layer_at_realign:end,:), align_shift, 2);
    ind_M(layer_at_realign:end,:) = circshift(ind_M(layer_at_realign:end,:), align_shift, 2);
    skl_M(layer_at_realign:end,:) = circshift(skl_M(layer_at_realign:end,:), align_shift, 2);
    
    
    % [debug] keep a record of trucated but not v-direction-stretched
    pxl_M_truncated = zeros(size(pxl_M));
    ind_M_truncated = zeros(size(ind_M));
    skl_M_truncated = zeros(size(skl_M));
    
    % This is to make each colomn stop at the first skeleton point it encounters, and then scale
    [nr,nc] = size(ind_M);
    for ic = 1:nc
        ind_r = find(skl_M(:,ic), 1, 'first');
        if isempty(ind_r)
            ind_r = size(skl_M,1);
        end
        
        % [debug] keep a record of truncated but not v-direction-stretched
        pxl_M_truncated(1:ind_r,ic) = pxl_M(1:ind_r,ic);
        ind_M_truncated(1:ind_r,ic) = ind_M(1:ind_r,ic);
        skl_M_truncated(1:ind_r,ic) = skl_M(1:ind_r,ic);
        
        % I am going to try two methods to interp, but I don't know which one is better.  This needs double check later.
        interpMethod = 1;
        if interpMethod == 1
            % in case skeleton is met at first row
            if ind_r==1
                ind_r=2;
            end
            pxl_M(:,ic) = interp1(1:ind_r, pxl_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
            ind_M(:,ic) = interp1(1:ind_r, ind_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
            skl_M(:,ic) = interp1(1:ind_r, skl_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
        elseif interpMethod == 2    % try to make only the last row as skeleton, not too much effective.
            if ind_r<=2
                pxl_M(1:end-1,ic) = interp1(1:2, pxl_M([1;1],ic), linspace(1,2,nr-1), 'nearest');
                ind_M(1:end-1,ic) = interp1(1:2, ind_M([1;1],ic), linspace(1,2,nr-1), 'nearest');
                skl_M(1:end-1,ic) = interp1(1:2, skl_M([1;1],ic), linspace(1,2,nr-1), 'nearest');
            else
                pxl_M(1:end-1,ic) = interp1(1:ind_r-1, pxl_M(1:ind_r-1,ic), linspace(1,ind_r-1,nr-1), 'nearest');
                ind_M(1:end-1,ic) = interp1(1:ind_r-1, ind_M(1:ind_r-1,ic), linspace(1,ind_r-1,nr-1), 'nearest');
                skl_M(1:end-1,ic) = interp1(1:ind_r-1, skl_M(1:ind_r-1,ic), linspace(1,ind_r-1,nr-1), 'nearest');
            end
            pxl_M(end,ic) = pxl_M(ind_r,ic);
            ind_M(end,ic) = ind_M(ind_r,ic);
            skl_M(end,ic) = ind_M(ind_r,ic);
        end
    end
else
    pxl_M = pxl_cell;
    ind_M = ind_cell;
    skl_M = skl_cell;
end

end

