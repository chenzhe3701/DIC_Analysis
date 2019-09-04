

% For each grain, summarize a [b#,g#] pair for distMap purpose, for each twin activity.     
% Also, generate one data point from the mantle, contribution to one of the previous defined regions.  

iE = 3;
eMap = calculate_effective_strain(strainFile{iE-1}.exx, strainFile{iE-1}.exy, strainFile{iE-1}.eyy);

% Find triple points ----------------------------------------------------
[boundary,~, ~, tripleTF, ~, indTriple, tripleIDs] = find_one_boundary_from_ID_matrix(ID);
xTriple = X(indTriple);
yTriple = Y(indTriple);

% ------------------------------------------------------------------------ 

% If exist data file, just load.  Otherwise, rerun analysis.  
try
    load(['temp_results/201908_quants_wrt_gb_',num2str(iE),'.mat']);
catch
    % summary of [bNum, gNum] pairs list, for different categories of activities
    % (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,
    % (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
    % (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth
    bg_not_involved = [];
    bg_slip_twin_a = [];
    bg_slip_twin_b =[];
    bg_co_found = [];
    bg_twin_twin_a = [];
    bg_twin_twin_b = [];
    bg_slip_growth_a = [];
    bg_slip_growth_b = [];
    bg_co_growth = [];
    
    % Summarize quants
    qs_not_involved = [];
    qs_slip_twin_a = [];
    qs_slip_twin_b = [];
    qs_co_found = [];
    qs_twin_twin_a = [];
    qs_twin_twin_b = [];
    qs_slip_growth_a = [];
    qs_slip_growth_b = [];
    qs_co_growth = [];
    
    NN = 250;
    x_dist = (1:NN)';
    qs_target = [0.0014, 0.0227, 0.1587, 0.5, 0.8414, 0.9772, 0.9987];
    
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        ind = find(gID==ID_current);
        nNeighbors = gNNeighbors(ind);
        ID_neighbors = gNeighbors(ind, 1:nNeighbors);
        
        ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        % Make it one data point wider on each side
        indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
        indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
        indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
        indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
        
        ID_local = ID(indR_min:indR_max, indC_min:indC_max);
        eMap_local = eMap(indR_min:indR_max, indC_min:indC_max);
        
        for iNb = 1:nNeighbors
            ID_neighbor = ID_neighbors(iNb);
            
            if ID_current < ID_neighbor
                gb = 10000*ID_neighbor + ID_current;
            else
                gb = 10000*ID_current + ID_neighbor;
            end
            
            
            
            % -------------------->       
            % Note, for the new analysis, we do not only want to find out what are the twin systems at this gb [iE, iTwin], we also want to find out [iE, iTwin, dist_to_triple]   
            % find the related triple points, which are surrounded by both ID_current and ID_neighbor
            ind_of_indTriple = sum(ismember(tripleIDs,ID_current) + ismember(tripleIDs,ID_neighbor),2)==2;
            x_triple = X(indTriple(ind_of_indTriple));
            y_triple = Y(indTriple(ind_of_indTriple));
            
            % A function to return activities in this grain, at this boundary, up to this iE (twin# intersect this gb at iEs)
            [iE_iTwin_dToTriple_list, valid_grain_a] = find_activity_at_boundary_from_struCell(struCell, iE, ID_current, gb);
            [iE_iTwin_dToTriple_list_Nb,valid_grain_b] = find_activity_at_boundary_from_struCell(struCell, iE, ID_neighbor, gb);
            
            if (valid_grain_a&&valid_grain_b)
                iE_g = min(iE_iTwin_dToTriple_list(:,1));
                iE_n = min(iE_iTwin_dToTriple_list_Nb(:,1));
                
                mm = 0;
                MM = 150;
                distMap_local = distance_from_boundary_in_grain(ID_local, [gb,ID_current]);
                ind = (distMap_local>mm)&(distMap_local<MM);
                
                if (iE_g>iE)
                    if (iE_n>iE)
                        % not involved
                        bg_not_involved = [bg_not_involved; gb, ID_current];
                        qs_not_involved = [qs_not_involved; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n==iE)
                        % slip induced twin, slip side
                        bg_slip_twin_a = [bg_slip_twin_a; gb, ID_current];
                        qs_slip_twin_a = [qs_slip_twin_a; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n<iE)
                        % slip induced twin growth, slip side
                        bg_slip_growth_a = [bg_slip_growth_a; gb, ID_current];
                        qs_slip_growth_a = [qs_slip_growth_a; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    end
                elseif (iE_g==iE)
                    % just twinned this grain
                    if (iE_n>iE)
                        % slip induced twin, twin side
                        bg_slip_twin_b = [bg_slip_twin_b; gb, ID_current];
                        qs_slip_twin_b = [qs_slip_twin_b; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n==iE)
                        % co-found twin
                        bg_co_found = [bg_co_found; gb, ID_current];
                        qs_co_found = [qs_co_found; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n<iE)
                        % twin induced twin, new twin side
                        bg_twin_twin_b = [bg_twin_twin_b; gb, ID_current];
                        qs_twin_twin_b = [qs_twin_twin_b; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    end
                elseif (iE_g<iE)
                    % twinned this grain
                    if (iE_n>iE)
                        % slip induced twin growth, twin side
                        bg_slip_growth_b = [bg_slip_growth_b; gb, ID_current];
                        qs_slip_growth_b = [qs_slip_growth_b; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n==iE)
                        % twin induced twin, old twin side
                        bg_twin_twin_a = [bg_twin_twin_a; gb, ID_current];
                        qs_twin_twin_a = [qs_twin_twin_a; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    elseif (iE_n<iE)
                        % co-growth
                        bg_co_growth = [bg_co_growth; gb, ID_current];
                        qs_co_growth = [qs_co_growth; quantile(eMap_local(ind), qs_target), nanmean(eMap_local(ind))];
                    end
                end
            end
            
                        
        end
        
    end
    %% Generate distance maps and calculate statistics for data points that belong to each category 
    distMap_not_involved = distance_from_boundary_in_grain(ID, bg_not_involved);
    distMap_slip_twin_a = distance_from_boundary_in_grain(ID, bg_slip_twin_a);
    distMap_slip_twin_b = distance_from_boundary_in_grain(ID, bg_slip_twin_b);
    distMap_co_found = distance_from_boundary_in_grain(ID, bg_co_found);
    distMap_twin_twin_a = distance_from_boundary_in_grain(ID, bg_twin_twin_a);
    distMap_twin_twin_b = distance_from_boundary_in_grain(ID, bg_twin_twin_b);
    distMap_slip_growth_a = distance_from_boundary_in_grain(ID, bg_slip_growth_a);
    distMap_slip_growth_b = distance_from_boundary_in_grain(ID, bg_slip_growth_b);
    distMap_co_growth = distance_from_boundary_in_grain(ID, bg_co_growth);
    
    q_not_involved = (cell2mat(...
        arrayfun(@(x) [quantile(eMap(distMap_not_involved(:)==x), qs_target), nanmean(eMap(distMap_not_involved(:)==x))],...
        x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 1');
    q_slip_twin_a = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_twin_a(:)==x), qs_target), nanmean(eMap(distMap_slip_twin_a(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 2');
    q_slip_twin_b = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_twin_b(:)==x), qs_target), nanmean(eMap(distMap_slip_twin_b(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 3');
    q_co_found = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_co_found(:)==x), qs_target), nanmean(eMap(distMap_co_found(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 4');
    q_twin_twin_a = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_twin_twin_a(:)==x), qs_target), nanmean(eMap(distMap_twin_twin_a(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 5');
    q_twin_twin_b = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_twin_twin_b(:)==x), qs_target), nanmean(eMap(distMap_twin_twin_b(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 6');
    q_slip_growth_a = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_growth_a(:)==x), qs_target), nanmean(eMap(distMap_slip_growth_a(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 7');
    q_slip_growth_b = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_slip_growth_b(:)==x), qs_target), nanmean(eMap(distMap_slip_growth_b(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 8');
    q_co_growth = (cell2mat(arrayfun(@(x) [quantile(eMap(distMap_co_growth(:)==x), qs_target), nanmean(eMap(distMap_co_growth(:)==x))], x_dist, 'uniformoutput',0)))';
    disp('sets of quants calculated: 9');
    %% And save results
    timeStr = datestr(now,'yyyymmdd_HHMM');
    save([timeStr,'_quants_wrt_gb_',num2str(iE),'.mat'], 'x_dist',...
        'bg_not_involved','bg_co_found','bg_co_growth','bg_slip_twin_a','bg_slip_twin_b','bg_twin_twin_a','bg_twin_twin_b','bg_slip_growth_a','bg_slip_growth_b',...
        'q_not_involved','q_co_found','q_co_growth','q_slip_twin_a','q_slip_twin_b','q_twin_twin_a','q_twin_twin_b','q_slip_growth_a','q_slip_growth_b',...
        'distMap_not_involved','distMap_co_found','distMap_co_growth','distMap_slip_twin_a','distMap_slip_twin_b','distMap_twin_twin_a','distMap_twin_twin_b','distMap_slip_growth_a','distMap_slip_growth_b');

end

%% can plot distance map to illustrate
if 0
    myplot(distMap_not_involved,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_slip_twin_a,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_slip_twin_b,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_co_found,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_twin_twin_a,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_twin_twin_b,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_slip_growth_a,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_slip_growth_b,boundaryTFB + logical(trueTwinMapCell{iE}));
    myplot(distMap_co_growth,boundaryTFB + logical(trueTwinMapCell{iE}));    
end

%% plot the (e.g.,) median/mean of strain at distances=1,2,3,... to grain boundary, for different categories
ith = 4;
figure; hold on;
nPts = length(x_dist);    % or 150
% ndp_to_plot = 150;
plot(x_dist(1:nPts),q_not_involved(ith,1:nPts),'-','color',[0 0.6 0],'LineWidth',2);
plot(x_dist(1:nPts),q_slip_twin_a(ith,1:nPts),'-b','LineWidth',2);
plot(x_dist(1:nPts),q_slip_twin_b(ith,1:nPts),'--b','LineWidth',2);

plot(x_dist(1:nPts),q_co_found(ith,1:nPts),'-k','LineWidth',2);
plot(x_dist(1:nPts),q_twin_twin_a(ith,1:nPts),'-r','LineWidth',2);
plot(x_dist(1:nPts),q_twin_twin_b(ith,1:nPts),'--r','LineWidth',2);

plot(x_dist(1:nPts),q_slip_growth_a(ith,1:nPts),'-m','LineWidth',2);
plot(x_dist(1:nPts),q_slip_growth_b(ith,1:nPts),'--m','LineWidth',2);
plot(x_dist(1:nPts),q_co_growth(ith,1:nPts),'-c','LineWidth',2);

legend({'non-involved','slip-twin slip-side','slip-twin new-twin-side', 'co-found','twin-twin old-twin-side','twin-twin new-twin-side', 'slip-growth slip-side','slip-growth twin-side','co-growth'});
set(gca,'fontsize',16);
xlabel('Distance to grain boundary');
ylabel('Effective strain');
title(['iE=',num2str(iE)],'fontweight','normal');
% set(gca,'ylim',[0.015 0.04]);
%% --> maybe check, find out the points that represent the median.  See where it is located.

%% boxplot to summarize the strain distribution in a mantle of [0,150] data point size, for different types of regions.
% (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,
% (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
% (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth 
close all;
mm = 0;
MM = 150;
v = [];
g = [];
ind = (distMap_not_involved>mm)&(distMap_not_involved<MM);
v = [v; eMap(ind)];
g = [g; 1*ones(sum(ind(:)),1)];
ind = (distMap_slip_twin_a>mm)&(distMap_slip_twin_a<MM);
v = [v; eMap(ind)];
g = [g; 2*ones(sum(ind(:)),1)];
ind = (distMap_slip_twin_b>mm)&(distMap_slip_twin_b<MM);
v = [v; eMap(ind)];
g = [g; 3*ones(sum(ind(:)),1)];

ind = (distMap_co_found>mm)&(distMap_co_found<MM);
v = [v; eMap(ind)];
g = [g; 4*ones(sum(ind(:)),1)];
ind = (distMap_twin_twin_a>mm)&(distMap_twin_twin_a<MM);
v = [v; eMap(ind)];
g = [g; 5*ones(sum(ind(:)),1)];
ind = (distMap_twin_twin_b>mm)&(distMap_twin_twin_b<MM);
v = [v; eMap(ind)];
g = [g; 6*ones(sum(ind(:)),1)];

ind = (distMap_slip_growth_a>mm)&(distMap_slip_growth_a<MM);
v = [v; eMap(ind)];
g = [g; 7*ones(sum(ind(:)),1)];
ind = (distMap_slip_growth_b>mm)&(distMap_slip_growth_b<MM);
v = [v; eMap(ind)];
g = [g; 8*ones(sum(ind(:)),1)];
ind = (distMap_co_growth>mm)&(distMap_co_growth<MM);
v = [v; eMap(ind)];
g = [g; 9*ones(sum(ind(:)),1)];

figure;
boxplot(v,g,'Notch','on','symbol','');
set(gca,'ylim',[-0.002, 0.03]);


%% For each grain, generate one data point from the mantle, contribution to one of the previous defined regions.  Could combine, but maybe better to run each section for each method.
% (1) not involved, (2) slip-induced-twin, slip side, (3) slip-induced-twin, twin side,  
% (4) co-found, (5) twin-induced-twin, old twin side, (6) twin-induced-twin, new twin side,
% (7) slip-induced-growth, slip side, (8) slip-induced-growth, twin side, (9) co-growth 
ith = 8;
vv = [];
gg = [];
vv = [vv;qs_not_involved(:,ith)];
gg = [gg; 1 * ones(size(qs_not_involved,1),1)];
vv = [vv;qs_slip_twin_a(:,ith)];
gg = [gg; 2 * ones(size(qs_slip_twin_a,1),1)];
vv = [vv;qs_slip_twin_b(:,ith)];
gg = [gg; 3 * ones(size(qs_slip_twin_b,1),1)];
vv = [vv;qs_co_found(:,ith)];
gg = [gg; 4 * ones(size(qs_co_found,1),1)];
vv = [vv;qs_twin_twin_a(:,ith)];
gg = [gg; 5 * ones(size(qs_twin_twin_a,1),1)];
vv = [vv;qs_twin_twin_b(:,ith)];
gg = [gg; 6 * ones(size(qs_twin_twin_b,1),1)];
vv = [vv;qs_slip_growth_a(:,ith)];
gg = [gg; 7 * ones(size(qs_slip_growth_a,1),1)];
vv = [vv;qs_slip_growth_b(:,ith)];
gg = [gg; 8 * ones(size(qs_slip_growth_b,1),1)];
vv = [vv;qs_co_growth(:,ith)];
gg = [gg; 9 * ones(size(qs_co_growth,1),1)];

figure;
boxplot(vv,gg,'Notch','on');
set(gca,'ylim',[0.004, 0.02]);
set(gca,'xticklabel',{'not involved','slip-twin slip side','slip-twin twin side',...
    'co-found','twin-twin, old twins side','twin-twin, new twin side',...
    'slip-growth slip side','slip-growth, twin side','co-growth'},'xticklabelrotation',45);


%%





