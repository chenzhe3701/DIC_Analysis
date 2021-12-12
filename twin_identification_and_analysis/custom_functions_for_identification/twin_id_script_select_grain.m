
% This is the script to modify twin identification. For different strain levels.
for ii = iE_start:iE_stop
    try
        close(ff{ii});
    catch
    end
end

% Select the grain of interest.
ids = find_ID_on_map(X,Y,ID,gcf,gca);
ID_current = ids(1);

% (1) Calculate theoretical trace direction.
ind_euler = find(gID==ID_current);
euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
if (1==eulerAligned)
    % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
    [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
else
    % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
    [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
end
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
mult_factor = ones(size(traceDir));
mult_factor(traceDir<0) = 1;
mult_factor(traceDir>=0) = -1;
traceND = traceDir + 90*mult_factor;    % convert traceDir to traceND
traceSF = abs_schmid_factor(nss+1:nss+ntwin,2);

disp(table(traceND,traceSF,'VariableNames',{'traceND','traceSF'}));

% (2) Here, we want to do analysis on clustered maps, to decide the active slip/twin system.
% For examples in twinning analysis, we previously performed cluster analysis and 'identified'/'confirmed' twin clusters.
% More generally, we might need to first do a rough clustering based on strain map, then perform trace analysis, to decide which are the exist slip/twin systems.

ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
indC_min = find(sum(ind_local, 1), 1, 'first');
indC_max = find(sum(ind_local, 1), 1, 'last');
indR_min = find(sum(ind_local, 2), 1, 'first');
indR_max = find(sum(ind_local, 2), 1, 'last');

ID_local = ID(indR_min:indR_max, indC_min:indC_max);

boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
x_local = X(indR_min:indR_max, indC_min:indC_max);
y_local = Y(indR_min:indR_max, indC_min:indC_max);
%%
for iE = iE_start:iE_stop
try
    close(ff{iE});
catch
end
end
for iE = iE_start:iE_stop
    ff{iE} = figure;
    eMapLocal = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
    tMapLocal = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    cMapLocal = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    
    eMapLocal(ID_local~=ID_current) = nan;
    tMapLocal(ID_local~=ID_current) = nan;
    cMapLocal(ID_local~=ID_current) = nan;
    
    eMapLocal(1==boundaryTF_local) = -inf;
    tMapLocal(1==boundaryTF_local) = -inf;
    cMapLocal(1==boundaryTF_local) = -inf;
        
    subplot(2,2,1);
    alpha = ones(size(eMapLocal));
    alpha(isnan(eMapLocal)) = 0;
    clim = quantile(eMapLocal(:),[0.005, 0.995]);
    imagesc([x_local(1),x_local(end)],[y_local(1),y_local(end)],eMapLocal,'alphadata',alpha);
    colorbar;
    try
        caxis(gca, [clim(1), clim(2)]);
    end
    axis equal;
    
    
    subplot(2,2,2);
    alpha = ones(size(tMapLocal));
    alpha(isnan(tMapLocal)) = 0;
    imagesc([x_local(1),x_local(end)],[y_local(1),y_local(end)],tMapLocal,'alphadata',alpha);
    colorbar;
    caxis([18,24]); 
    axis equal;
    
    
    subplot(2,2,3);
    alpha = ones(size(cMapLocal));
    alpha(isnan(cMapLocal)) = 0;
    imagesc([x_local(1),x_local(end)],[y_local(1),y_local(end)],cMapLocal,'alphadata',alpha);
    colorbar;
%     caxis([0, 5]);
    axis equal;
    
end








