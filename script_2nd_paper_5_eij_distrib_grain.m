% chenzhe, 2018-04-07, add note
% move code from code '..._twin_tracking_example' to here as individual code  
%
% continue on 2018-03-21, examine pointwise distance to twin strain
% try using K=2-5 clusters, plot 
% (1) pointwise distance of strain components to highest-SFed twin center. 
% This could be used to make a contour plot of equi-distance, but it turned
% out hard to see, so just use this pointwise distance map.
% (2) The histogram of this distance. The boundary of twin vs non-twin can
% be guessed from the position between peaks. -- Important thing is, the
% boundary position does not seem to be constant.

% First, iE of interest ----------------------------- !!! -------------------------------------

[fname, grainDataPath] = uigetfile('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\g_1144_new_20180321\WE43_T6_C1_s_all_grain_1144_local_map.mat', 'select grain data');

load(fullfile(grainDataPath,fname),'data');
clear exx_distribution exy_distribution eyy_distribution

for iE = 2:5
    s = data(iE);
    
    indR_min = s.indR_min;
    indR_max = s.indR_max;
    indC_min = s.indC_min;
    indC_max = s.indC_max;
    boundaryTF_local = s.boundaryTF_local;
    x_local = s.x_local;
    y_local = s.y_local;
    ID_current = s.ID_current;
    ID_local = s.ID_local;
    exx_local = s.exx_local;
    exy_local = s.exy_local;
    eyy_local = s.eyy_local;
    sigma_local = s.sigma_local;
    
    % summarize the distribution of exx data at different iE level
    % figure;histogram(exx_local(:));
    % figure;histogram(exy_local(:));
    % figure;histogram(eyy_local(:));
    edges_xx = linspace(-0.2,0.05,100);
    edges_xy = linspace(-0.04,0.06,100);
    edges_yy = linspace(-0.04,0.06,100);
    exx_distribution{iE} = histcounts(exx_local(:),edges_xx);
    exy_distribution{iE} = histcounts(exy_local(:),edges_xx);
    eyy_distribution{iE} = histcounts(eyy_local(:),edges_xx);
    
 end

%% After the previous section, plot distribution of exx, exy, and eyy
close all;
colors = lines(7);

figure; hold on;
plot(edges_xx(1:end-1),exx_distribution{3},'-','color',colors(1,:),'linewidth',1.5);
plot(edges_xx(1:end-1),exx_distribution{4},'-','color',colors(2,:),'linewidth',1.5);
plot(edges_xx(1:end-1),exx_distribution{5},'-','color',colors(4,:),'linewidth',1.5);

legend({'Global Uniaxial Strain -1.2%','Global Uniaxial Strain -2.1%','Global Uniaxial Strain -3.7%'},'location','northwest');
xlabel('\epsilon_x_x value, mm/mm');
ylabel('Counts');
set(gca,'fontsize',16);
print(fullfile(grainDataPath,'exx_distrib_curve.tif'),'-dtiff');


figure; hold on;
plot(edges_xy(1:end-1),exy_distribution{3},'-','color',colors(1,:),'linewidth',1.5);
plot(edges_xy(1:end-1),exy_distribution{4},'-','color',colors(2,:),'linewidth',1.5);
plot(edges_xy(1:end-1),exy_distribution{5},'-','color',colors(4,:),'linewidth',1.5);

legend({'Global Uniaxial Strain -1.2%','Global Uniaxial Strain -2.1%','Global Uniaxial Strain -3.7%'},'location','northwest');
xlabel('\epsilon_x_y value, mm/mm');
ylabel('Counts');
set(gca,'fontsize',16)
print(fullfile(grainDataPath,'exy_distrib_curve.tif'),'-dtiff');

figure; hold on;
plot(edges_yy(1:end-1),eyy_distribution{3},'-','color',colors(1,:),'linewidth',1.5);
plot(edges_yy(1:end-1),eyy_distribution{4},'-','color',colors(2,:),'linewidth',1.5);
plot(edges_yy(1:end-1),eyy_distribution{5},'-','color',colors(4,:),'linewidth',1.5);

legend({'Global Uniaxial Strain -1.2%','Global Uniaxial Strain -2.1%','Global Uniaxial Strain -3.7%'},'location','northwest');
xlabel('\epsilon_y_y value, mm/mm');
ylabel('Counts');
set(gca,'fontsize',16)
print(fullfile(grainDataPath,'eyy_distrib_curve.tif'),'-dtiff');




