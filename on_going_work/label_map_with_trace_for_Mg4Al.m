
% Need to have euler angles in base
% currently only for Mg.
%
% I can improve it later when I have time.  chenzhe, 2019-08-02

function [sfs,angles] = label_map_with_trace(X,Y,ID,ID_targets,ss_target,ha)

gID = evalin('base','gID;');
gPhi1 = evalin('base','gPhi1;');
gPhi = evalin('base','gPhi;');
gPhi2 = evalin('base','gPhi2;');
eulerAligned = evalin('base','eulerAligned;');
stressTensor = evalin('base','stressTensor;');
sampleMaterial = evalin('base','sampleMaterial;');

% get stats about neighbor and plot, such as m'
try
    ID_reduced = ID(1:20:end,1:20:end);
    X_reduced = X(1:20:end,1:20:end);
    Y_reduced = Y(1:20:end,1:20:end);
    unique(ID_reduced(:)) - unique(ID(:));
catch
    ID_reduced = ID(1:3:end,1:3:end);
    X_reduced = X(1:3:end,1:3:end);
    Y_reduced = Y(1:3:end,1:3:end);
end

for ii=1:length(ID_targets)
    ID_target = ID_targets(ii);
    ind = find(gID==ID_target);
    euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
    % Now I want to
    % (1) label basal slip trace
    if (1==eulerAligned)
        % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
        [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
    else
        % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
        [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
    end
    
    % [max_basal_sf, ind_ss] = max(abs_schmid_factor(ss_target(:),2));
    % angle = abs_schmid_factor(ind_ss,3);
    ssLabel = abs_schmid_factor(:,1);
    sfs = abs_schmid_factor(:,2);
    angles = abs_schmid_factor(:,3);
    
    inds = (ID_reduced==ID_target);
    x_this_nb = mean(X_reduced(inds));
    y_this_nb = mean(Y_reduced(inds));
    x_range = range(X_reduced(inds));
    y_range = range(Y_reduced(inds));
    
    try
        % plot from low SF to high SF for color purpose. --------
        for ii = 1:length(ss_target)
            is = ss_target(ii);
            
            angle = angles(is);
            dx = x_range/4;
            dy = x_range/4 * tand(angle);
            if abs(dy) > y_range/4
                factor = y_range/4 / dy;
                dy = dy * factor;
                dx = dx * factor;
            end
            
            if is<=3
                plot3(ha, [x_this_nb - dx, x_this_nb + dx], [y_this_nb - dy, y_this_nb + dy], [20,20],'linewidth',1,'color', 'r');
            elseif is<=6
                plot3(ha, [x_this_nb - dx, x_this_nb + dx], [y_this_nb - dy, y_this_nb + dy], [20,20],'linewidth',1,'color', 'k');
            elseif is<=12
            elseif is<=18
            else
                plot3(ha, [x_this_nb - dx, x_this_nb + dx], [y_this_nb - dy, y_this_nb + dy], [20,20],'linewidth',1,'color', 'g');
            end
            %         text(x_this_nb, y_this_nb, 50, num2str(ID_target));
            %     text(x_this_nb, y_this_nb, 50, num2str(sfs(ii),3));
        end
    catch
        disp(['error in finding ID: ',num2str(ID_target)]);
    end
    
end

