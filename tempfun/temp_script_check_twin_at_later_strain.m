
% This is just check the twin label at different strain level.
% Use as a helper for 'script_manual_label_twin_gb_intersection.m'

iE 

ID_current = struCell{iE}(iS).gID;
ind = find(gID==ID_current);

euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
if (1==eulerAligned)
    % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
    [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
else
    % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
    [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
end
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
% mult_factor = ones(size(traceDir));
% mult_factor(traceDir<0) = 1;
% mult_factor(traceDir>=0) = -1;
% traceND = traceDir + 90*mult_factor;    % convert traceDir to traceND


nNeighbors = gNNeighbors(ind);
ID_neighbors = gNeighbors(ind, 1:nNeighbors);

ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);

% Make it one data point wider on each side
indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);

ID_local = ID(indR_min:indR_max, indC_min:indC_max);
X_local = X(indR_min:indR_max, indC_min:indC_max);
Y_local = Y(indR_min:indR_max, indC_min:indC_max);
uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
trueTwinMapLocal = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
e_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);

[handleFig0,a1,~] = myplot(X_local, Y_local, e_local, grow_boundary(boundaryTF_local));
disableDefaultInteractivity(a1);
[handleFig,aa,~] = myplot(X_local, Y_local, trueTwinMapLocal, grow_boundary(boundaryTF_local));
caxis([18 24]);
label_map_with_ID(X_local, Y_local, ID_local, handleFig, ID_current);
disableDefaultInteractivity(aa);