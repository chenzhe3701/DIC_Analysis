% [pxl_M, ind_M, skl_M] = onion_thin(A, ind_start_input)
% A = input image
% ind_start_point is the index to start traversal.
% chenzhe, 2018-05-19



function [pxl_M, ind_M, skl_M] = onion_thin(A, ind_start_input)

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

pxl_val_cell = [];
ind_cell = [];
skl_cell = [];

% start at
ind = ind_start_input;
reAligned = false;

%% Enter level-1 loop, for each layer
all_skeleton = 0;
iLoop = 1;
while (~all_skeleton)&&(iLoop<=inf)
    
    pxl_val_list = [];
    ind_list = [];
    skl_list = [];
    
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
        pxl_val_list = [pxl_val_list, A_raw(ind)];
        ind_list = [ind_list, ind]; 
        skl_list = [skl_list, skeleton(ind)];
        
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
        disp(['Finished traversing at Loop: ', num2str(iLoop)]);
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
        else
            % if only have one point in this layer, make it 2 points to interp 
            if length(ind_list)==1
                pxl_val_list = [pxl_val_list,pxl_val_list];
                ind_list = [ind_list,ind_list];
                skl_list = [skl_list,skl_list];
            end
            pxl_val_list = interp1(1:length(pxl_val_list), pxl_val_list, linspace(1,length(pxl_val_list),list_length),'nearest');
            ind_list = interp1(1:length(ind_list), ind_list, linspace(1,length(ind_list),list_length),'nearest');
            skl_list = interp1(1:length(skl_list), skl_list, linspace(1,length(skl_list),list_length),'nearest');
        end
        
        % find starting point of next layer
        ind = nb_ind_list(find(A(nb_ind_list)==1, 1, 'first'));
        
        % try to see if can re-align
        align_shift = find(skl_list, 1, 'first')-1;
        
        % It can be re-aligned
        if ~isempty(align_shift)
            reAligned = true;
            pxl_val_list = circshift(pxl_val_list, -align_shift, 2);
            ind_list = circshift(ind_list, -align_shift, 2);
            skl_list = circshift(skl_list, -align_shift, 2);
            layer_at_realign = iLoop;
            ind = ind_list(1);
            ind_previous = ind_list(end);
        end
        % copy the list into cell. If realigned, then copy the realigned version, and record info: [align_shift, layer_at_realign]     
        pxl_val_cell{iLoop} = pxl_val_list;
        ind_cell{iLoop} = ind_list;
        skl_cell{iLoop} = skl_list;     
    else
        % align by considering anchor points
        % Use last loops aligned ind_list as template, find all the anchor positions in index in this list, make sure the last one is considered as anchor  
        % Only interp between different GROUPs of anchor points.
        % --> If two anchor points are in the same group in step-i, they will remain in the same group in step i+1.
        % It is possible that some anchor points in a group in step i can show twice in step i+1.
        % But if we interp using the largest range, because anchors deccendants will be eliminated, so it is still OK to just interp.  
        % I think by finding the 'last' of p2, it can be achieved that two anchors are from different 'groups'  
        % example-1: 
        % layer-1: [14, 15, 16, ... ...]
        % layer-2: [14, 16, 14, ...]
        % example-2:
        % layer-1: [226, 225, ..., ..., ..., ..., ..., 715]
        % layer-2: [226, 225, 715, ..., 265, 715, 225]
        anchor_ind = ind_cell{iLoop-1};
        anchor_pos = find(skl_cell{iLoop-1});
        if anchor_pos(end)~=list_length
            anchor_pos = [anchor_pos,list_length];
        end
        % interp at different stepSize at each anchor interval
        xp = [3.141592653589793];   % request point, intialize in order to have an 'xp(end)', -- the value does not matter.
        for ii=1:length(anchor_pos)-1
           p1 = find(ind_list==anchor_ind(anchor_pos(ii)), 1, 'first');
           p2 = find(ind_list==anchor_ind(anchor_pos(ii+1)), 1, 'last');
           % if the last point in ind_anchor is not a skeleton point, assign the position as anchor  
           if isempty(p2) && (ii==length(anchor_pos)-1)
              p2 = length(ind_list);
           end
           xp(end) = [];
           xp = [xp, linspace(p1,p2,anchor_pos(ii+1)-anchor_pos(ii)+1)];
        end
  
        % interp
        pxl_val_list = interp1(1:length(pxl_val_list), pxl_val_list, xp,'nearest');
        ind_list = interp1(1:length(ind_list), ind_list, xp,'nearest');
        skl_list = interp1(1:length(skl_list), skl_list, xp,'nearest');
        
        % copy the list into cell
        pxl_val_cell{iLoop} = pxl_val_list;
        ind_cell{iLoop} = ind_list;
        skl_cell{iLoop} = skl_list;    
    end
    
    % increment loop
    iLoop = iLoop + 1;
end

%% make into matrix
pxl_M = cell2mat(reshape(pxl_val_cell,[],1));
ind_M = cell2mat(reshape(ind_cell,[],1));
skl_M = cell2mat(reshape(skl_cell,[],1));

% shift back
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
    
    % [debug] keep a record of truncated but not v-direction-stretched
    pxl_M_truncated(1:ind_r,ic) = pxl_M(1:ind_r,ic);
    ind_M_truncated(1:ind_r,ic) = ind_M(1:ind_r,ic);
    skl_M_truncated(1:ind_r,ic) = skl_M(1:ind_r,ic);
    
    % in case skeleton is met at first row
    if ind_r==1
        ind_r=2;
    end
    
    pxl_M(:,ic) = interp1(1:ind_r, pxl_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
    ind_M(:,ic) = interp1(1:ind_r, ind_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
    skl_M(:,ic) = interp1(1:ind_r, skl_M(1:ind_r,ic), linspace(1,ind_r,nr), 'nearest');
end






