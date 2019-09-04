
% Need to have euler angles in base

function [sfs,angles] = local_fun_label_map_with_ID_and_trace(X,Y,ID,ID_target,ss_target,ha)

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

ind_Nb = find(gID==ID_target);
euler_Nb = [gPhi1(ind_Nb),gPhi(ind_Nb),gPhi2(ind_Nb)];
% Now I want to
% (1) label basal slip trace
if (1==eulerAligned)
    % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
    [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_Nb, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
else
    % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
    [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_Nb, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
end

% [max_basal_sf, ind_ss] = max(abs_schmid_factor(ss_target(:),2));
% angle = abs_schmid_factor(ind_ss,3);
ssLabel = sorted_schmid_factor(:,1);
sfs = sorted_schmid_factor(:,2);
angles = sorted_schmid_factor(:,3);

inds = (ID_reduced==ID_target);
x_this_nb = mean(X_reduced(inds));
y_this_nb = mean(Y_reduced(inds));
x_range = range(X_reduced(inds));
y_range = range(Y_reduced(inds));

% plot from low SF to high SF for color purpose. --------
for ii = length(sfs):-1:1
    if ismember(ssLabel(ii), ss_target)
        angle = angles(ii);
        dx = x_range/4;
        dy = x_range/4 * tand(angle);
        if abs(dy) > y_range/4
            factor = y_range/4 / dy;
            dy = dy * factor;
            dx = dx * factor;
        end
        if ismember(ssLabel(ii),1:3)
            plot3(ha, [x_this_nb - dx, x_this_nb + dx], [y_this_nb - dy, y_this_nb + dy], [20,20],'linewidth',2,'color', [1 0 0] * max(sfs(ii),0)*2);
        elseif ismember(ssLabel(ii),19:24)
            plot3(ha, [x_this_nb - dx, x_this_nb + dx], [y_this_nb - dy, y_this_nb + dy], [20,20],'linewidth',2,'color', [0 0 1] * max(sfs(ii),0)*2);
        end
        text(x_this_nb, y_this_nb, 50, num2str(ID_target));
        %     text(x_this_nb, y_this_nb, 50, num2str(sfs(ii),3));
    end
end


