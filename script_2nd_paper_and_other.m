% script to make figure for paper: 
% Can be used for MSEA paper.
% based on script for paper 1. But can be for other use
close all;
iE = 4;
strain_file = ['D:\WE43_T6_C1_insitu_compression\stitched_DIC\_',num2str(iE),'.mat'];
translation_file = 'D:\WE43_T6_C1_insitu_compression\stitched_img\translations_searched_vertical_stop_0.mat'
ws_file = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat'
twin_file = 'D:\p\m\DIC_Analysis\twinMaps.mat';

addChenFunction;
load(strain_file,'exx','u');
load(translation_file);
load(ws_file,'ID','boundaryTF','boundaryTFB','X','Y');
load(twin_file,'twinMapCell');
%% plot 
% plot exx map
[f1,a1,c1]=myplotm(exx,'x',X,'y',Y,'tf',boundaryTFB,'r',1);
clim = caxis

c = colorbar;
caxis([-0.12, 0.02])
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x, mm/mm');

%% Make double axis
script_make_double_axis()
title(c,'\epsilon_x_x');

%% enlarged

ID_target = 1144;   % This is the grain for MSEA
ID_target = 190;    % This is the grain for my research statement
ID_target = 639;    % This is the grain for 2020 paper
% ID_target = find_ID_on_map(X,Y,ID,gcf,gca);

ID_current=ID_target;
freeSpace = 10;

ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
indC_min = find(sum(ind_local, 1), 1, 'first') - freeSpace;
indC_max = find(sum(ind_local, 1), 1, 'last') + freeSpace;
indR_min = find(sum(ind_local, 2), 1, 'first') - freeSpace;
indR_max = find(sum(ind_local, 2), 1, 'last') + freeSpace;

exx_local = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
x_local = X(indR_min:indR_max, indC_min:indC_max);
y_local = Y(indR_min:indR_max, indC_min:indC_max);
ID_local = ID(indR_min:indR_max, indC_min:indC_max);
boundary_local = boundaryTFB(indR_min:indR_max, indC_min:indC_max);
    

% plot local exx
[f3,a3,c3]=myplotm(exx_local,'x',x_local,'y',y_local,'tf',boundary_local,'r',1);
c = colorbar
caxis([-0.12 0.02]);
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('');
title(c,'\epsilon_x_x','fontweight','normal');
axis off;

%%
twinMap_local = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
% plot twinMap_local
twinMap_local(twinMap_local<19) = nan;
[f4,a4,c4]=myplotm(twinMap_local,'x',x_local,'y',y_local,'tf',boundary_local,'r',1);
c = colorbar
caxis([18 25]);
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('');
title(c,'Twin System #','fontweight','normal');
axis off;


%% all the print functions
print(['exxMap_iE_',num2str(iE),'.tif'],'-dtiff');
print(['exxMapLocal_iE_',num2str(iE),'.tif'],'-dtiff');
print(['twinMapLocal_iE_',num2str(iE),'.tif'],'-dtiff');

