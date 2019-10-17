
% If want to twin, (1) cannot have high basal_SF (2) cannot have low twin_SF    

close all;
iE = 5;

% (1) max basal SF of grain
edges = 0:0.05:0.5;

ind1 = (T.iE==iE)&(T.twinnedTF==0);
N_nt = histcounts(T.basal_SF(ind1), edges);
ind2 = (T.iE==iE)&(T.twinnedTF==1);
N_t = histcounts(T.basal_SF(ind2), edges);
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
figure;disableDefaultInteractivity(gca);
bar(xpos, [N_nt(:),N_t(:)], 1, 'stacked');
xlabel('Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent grains twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Basal SF, iE = ',num2str(iE)],'fontweight','normal');

% (2) basal SF, T2, for all variants
edges = 0:0.05:0.5;

ind1 = (T2.iE==iE)&(T2.vActiveTF==0);
N_nt = histcounts(T2.basal_SF(ind1), edges);
ind2 = (T2.iE==iE)&(T2.vActiveTF==1);
N_t = histcounts(T2.basal_SF(ind2), edges);
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
figure;disableDefaultInteractivity(gca);
bar(xpos, [N_nt(:),N_t(:)], 1, 'stacked');
xlabel('Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Basal SF, iE = ',num2str(iE)],'fontweight','normal');

% (3) max twin SF of grain
edges = -0.5:0.05:0.5;

ind1 = (T.iE==iE)&(T.twinnedTF==0);
N_nt = histcounts(T.twin_SF(ind1), edges);
ind2 = (T.iE==iE)&(T.twinnedTF==1);
N_t = histcounts(T.twin_SF(ind2), edges);
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
figure;disableDefaultInteractivity(gca);
bar(xpos, [N_nt(:),N_t(:)], 1, 'stacked');
xlabel('Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Twin SF, iE = ',num2str(iE)],'fontweight','normal');

% (4) twin SF, T2, for all variants
edges = -0.5:0.05:0.5;

ind1 = (T2.iE==iE)&(T2.vActiveTF==0);
N_nt = histcounts(T2.variant_SF(ind1), edges);
ind2 = (T2.iE==iE)&(T2.vActiveTF==1);
N_t = histcounts(T2.variant_SF(ind2), edges);
d_int = (edges(3)-edges(2))/2;
xpos = [edges(2)-d_int, edges(2:end-1)+d_int];
xstr = [];
for ii=1:length(edges)-1
    xstr{ii} = [num2str(edges(ii)),'-',num2str(edges(ii+1))];
end
figure;disableDefaultInteractivity(gca);
bar(xpos, [N_nt(:),N_t(:)], 1, 'stacked');
xlabel('Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Twin SF, iE = ',num2str(iE)],'fontweight','normal');

% (5) (x-axis) twin SF, (y-axis) basal SF, for grains with (red) and without (black) twin  
figure; hold on;
ind = (T.iE==iE)&(T.twinnedTF==1);
plot(T.twin_SF(ind), T.basal_SF(ind),'.r','markersize',16);
ind = (T.iE==iE)&(T.twinnedTF==0);
plot(T.twin_SF(ind), T.basal_SF(ind),'.k');
xlabel('twin SF');
ylabel('basal SF');

%% 
save('temp_results\saved_data.mat','placeHolder','-v7.3')
%%
save('temp_results\saved_data.mat','edmat_2','-append','-v7.3')
%%
D = matfile('temp_results\saved_data.mat')
%% note, check, use -v7.3 file with outlier_removed=1.

% calculate edmat: strain distribution matrix
% [iE, ID_current, gb, nanmean of eEff @ d=1:250 data point distance to that gb.  

for iE = 2:5
    
    eMap = calculate_effective_strain(strainFile{iE-1}.exx, strainFile{iE-1}.exy, strainFile{iE-1}.eyy);
    edmat = [];
    
    iS = 1;
    warning('off','MATLAB:table:RowsAddedExistingVars');
    continueTF = true;
    dToTriple_th = 5;       % eliminate intersection whose distance to triple point is smaller than this value
    dToTriple_th_to_label = dToTriple_th;  % label if distance of intersection to triple point is smaller than this value
    
    hW = waitbar(0, ['iE=',num2str(iE),' analyze each grain']);
    hN = length(struCell{iE});
    while (continueTF)&&(iS<=length(struCell{iE}))
        waitbar(iS/hN, hW);
        
        close all;
        ID_current = struCell{iE}(iS).gID
        ind = find(gID==ID_current);
        
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
        
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);  % This is for effective strain
        
        activeTS = cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb);
        tSF = struCell{iE}(iS).tSF;
        
        % If grain of interest is twinned (here it means: got a twin-gb intersection labeled)
        if any(cellfun(@(x) ~isempty(x), struCell{iE}(iS).tGb)) % sum(struCell{iE}(iS).cTrueTwin(:))>0
            
            % [[[[For each neighbor]  get stats about neighbor and plot, such as m'
            for iNb = 1:nNeighbors
                ID_neighbor = ID_neighbors(iNb);
                iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
                if ~isempty(iS_neighbor)
                    activeTS_Nb = cellfun(@(x) ~isempty(x), struCell{iE}(iS_neighbor).tGb);
                    
                    % (1.1) Calculate this_uniqueGB number.
                    if ID_current > ID_neighbor
                        gb = ID_current * 10000 + ID_neighbor;
                    else
                        gb = ID_neighbor * 10000 + ID_current;
                    end
                    
                    % strain calculation in area of interest.
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                    edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];
                    
                end
                % end of ~isempty(iS_neighbor)
                
            end
            
        else
            % If grain of interest is not twinned, this means (iE_g > iE)
            
            for iNb = 1:nNeighbors
                ID_neighbor = ID_neighbors(iNb);
                
                iS_neighbor = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
                if ~isempty(iS_neighbor)
                    % (1.1) Calculate this_uniqueGB number.
                    if ID_current > ID_neighbor
                        gb = ID_current * 10000 + ID_neighbor;
                    else
                        gb = ID_neighbor * 10000 + ID_current;
                    end
                    
                    % strain calculation in area of interest.
                    distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                    edmat = [edmat; iE, ID_current, gb, arrayfun(@(x) nanmean(eMap_local(distMap_local==x)), 1:250)];
                else
                    disp(['neighbor not found: ',num2str(ID_neighbor)]);
                end
            end
        end
        
        % disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
        iS = iS + 1;
    end
    close(hW);
    warning on;
    
    switch iE
        case 2
            edmat_2 = edmat;
            save('temp_results\saved_data.mat','edmat_2','-append','-v7.3');
        case 3
            edmat_3 = edmat;
            save('temp_results\saved_data.mat','edmat_3','-append','-v7.3');
        case 4
            edmat_4 = edmat;
            save('temp_results\saved_data.mat','edmat_4','-append','-v7.3');
        case 5
            edmat_5 = edmat;
            save('temp_results\saved_data.mat','edmat_5','-append','-v7.3');
    end
    
end

%%  Plot curve.
iE = 2;
load(['D:\p\m\DIC_Analysis\temp_results\quants_wrt_gb_',num2str(iE),'.mat']);
load('D:\p\m\DIC_Analysis\temp_results\saved_data.mat', ['edmat_',num2str(iE)]);

expr = ['edmat_',num2str(iE),';'];
edmat = evalin('base',expr);

% nanmean
clear eline_not_involved eline_slip_twin_a eline_slip_twin_b eline_co_found eline_twin_twin_a eline_twin_twin_b eline_slip_growth_a eline_slip_growth_b eline_co_growth
legend_str = [];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_not_involved, 'rows');
eline_not_involved = nanmean(edmat(ind,4:end),1);
legend_str{1} = ['not-involved: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_a, 'rows');
eline_slip_twin_a = nanmean(edmat(ind,4:end),1);
legend_str{2} = ['slip-twin slip-side: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_b, 'rows');
eline_slip_twin_b = nanmean(edmat(ind,4:end),1);
legend_str{3} = ['slip-twin new-twin-side: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_found, 'rows');
eline_co_found = nanmean(edmat(ind,4:end),1);
legend_str{4} = ['co-found: ',num2str(sum(ind))];
try
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_a, 'rows');
    eline_twin_twin_a = nanmean(edmat(ind,4:end),1);
    legend_str{5} = ['twin-twin old-twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_b, 'rows');
    eline_twin_twin_b = nanmean(edmat(ind,4:end),1);
    legend_str{6} = ['twin-twin new-twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_a, 'rows');
    eline_slip_growth_a = nanmean(edmat(ind,4:end),1);
    legend_str{7} = ['slip-growth slip-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_b, 'rows');
    eline_slip_growth_b = nanmean(edmat(ind,4:end),1);
    legend_str{8} = ['slip-growth twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_growth, 'rows');
    eline_co_growth = nanmean(edmat(ind,4:end),1);
    legend_str{9} = ['co-growth: ',num2str(sum(ind))];
end
figure;disableDefaultInteractivity(gca); hold on;
nPts = length(x_dist);    % or 150
% ndp_to_plot = 150;
plot(x_dist(1:nPts),eline_not_involved(1:nPts),'-','color',[0 0.6 0],'LineWidth',3);
plot(x_dist(1:nPts),eline_slip_twin_a(1:nPts),'-b','LineWidth',3);
plot(x_dist(1:nPts),eline_slip_twin_b(1:nPts),'--b','LineWidth',2);

plot(x_dist(1:nPts),eline_co_found(1:nPts),'-k','LineWidth',2);
try
    plot(x_dist(1:nPts),eline_twin_twin_a(1:nPts),'-r','LineWidth',.5);
    plot(x_dist(1:nPts),eline_twin_twin_b(1:nPts),'--r','LineWidth',2);
    
    plot(x_dist(1:nPts),eline_slip_growth_a(1:nPts),'-m','LineWidth',3);
    plot(x_dist(1:nPts),eline_slip_growth_b(1:nPts),'--m','LineWidth',.5);
    plot(x_dist(1:nPts),eline_co_growth(1:nPts),'-c','LineWidth',.5);
end
% legend({'non-involved','slip-twin slip-side','slip-twin new-twin-side', 'co-found','twin-twin old-twin-side','twin-twin new-twin-side', 'slip-growth slip-side','slip-growth twin-side','co-growth'});
legend(legend_str);
set(gca,'fontsize',12);
xlabel('Distance to grain boundary');
ylabel('Mean of Effective strain');

title(['iE=',num2str(iE)],'fontweight','normal');

% median
clear eline_not_involved eline_slip_twin_a eline_slip_twin_b eline_co_found eline_twin_twin_a eline_twin_twin_b eline_slip_growth_a eline_slip_growth_b eline_co_growth

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_not_involved, 'rows');
eline_not_involved = median(edmat(ind,4:end),1,'omitnan');
legend_str{1} = ['not-involved: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_a, 'rows');
eline_slip_twin_a = median(edmat(ind,4:end),1,'omitnan');
legend_str{2} = ['slip-twin slip-side: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_twin_b, 'rows');
eline_slip_twin_b = median(edmat(ind,4:end),1,'omitnan');
legend_str{3} = ['slip-twin new-twin-side: ',num2str(sum(ind))];

[ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_found, 'rows');
eline_co_found = median(edmat(ind,4:end),1,'omitnan');
legend_str{4} = ['co-found: ',num2str(sum(ind))];
try
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_a, 'rows');
    eline_twin_twin_a = median(edmat(ind,4:end),1,'omitnan');
    legend_str{5} = ['twin-twin old-twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_twin_twin_b, 'rows');
    eline_twin_twin_b = median(edmat(ind,4:end),1,'omitnan');
    legend_str{6} = ['twin-twin new-twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_a, 'rows');
    eline_slip_growth_a = median(edmat(ind,4:end),1,'omitnan');
    legend_str{7} = ['slip-growth slip-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_slip_growth_b, 'rows');
    eline_slip_growth_b = median(edmat(ind,4:end),1,'omitnan');
    legend_str{8} = ['slip-growth twin-side: ',num2str(sum(ind))];
    
    [ind,~] = ismember([edmat(:,3),edmat(:,2)], bg_co_growth, 'rows');
    eline_co_growth = median(edmat(ind,4:end),1,'omitnan');
    legend_str{9} = ['co-growth: ',num2str(sum(ind))];
end

figure;disableDefaultInteractivity(gca); hold on;
nPts = length(x_dist);    % or 150
% ndp_to_plot = 150;
plot(x_dist(1:nPts),eline_not_involved(1:nPts),'-','color',[0 0.6 0],'LineWidth',3);
plot(x_dist(1:nPts),eline_slip_twin_a(1:nPts),'-b','LineWidth',3);
plot(x_dist(1:nPts),eline_slip_twin_b(1:nPts),'--b','LineWidth',2);

plot(x_dist(1:nPts),eline_co_found(1:nPts),'-k','LineWidth',2);
try
    plot(x_dist(1:nPts),eline_twin_twin_a(1:nPts),'-r','LineWidth',.5);
    plot(x_dist(1:nPts),eline_twin_twin_b(1:nPts),'--r','LineWidth',2);
    
    plot(x_dist(1:nPts),eline_slip_growth_a(1:nPts),'-m','LineWidth',3);
    plot(x_dist(1:nPts),eline_slip_growth_b(1:nPts),'--m','LineWidth',.5);
    plot(x_dist(1:nPts),eline_co_growth(1:nPts),'-c','LineWidth',.5);
end
% legend({'non-involved','slip-twin slip-side','slip-twin new-twin-side', 'co-found','twin-twin old-twin-side','twin-twin new-twin-side', 'slip-growth slip-side','slip-growth twin-side','co-growth'});
legend(legend_str);
set(gca,'fontsize',12);
xlabel('Distance to grain boundary');
ylabel('Medium of Effective strain');

title(['iE=',num2str(iE)],'fontweight','normal');









