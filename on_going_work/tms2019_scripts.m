%% Use together with script_twinVariantAnalysisProcedure.m

%% print exx map for grain 694 at ie=4, 
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','')
title('');
print('figFolder\g694_ie4_exx.tiff','-dtiff')

%% print cluster number map raw
colormap(parula(3))
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
caxis([-0.5 2.5]);
c=colorbar;
set(c,'limits',[0.5 2.5],'Ticks',[1 2])
title('');
print('figFolder\g694_ie4_cNumMap_raw.tiff','-dtiff');

%% print cluster number map cleaned
colormap(parula(3))
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
caxis([-0.5 2.5]);
c=colorbar;
set(c,'limits',[0.5 2.5],'Ticks',[1 2])
title('');
print('figFolder\g694_ie4_cNumMap_cleaned.tiff','-dtiff');

%% print skeleton map
colorbar off;
title('');
print('figFolder\g694_ie4_skeleton.tiff','-dtiff');

%% print the Hough transform map
set(gca,'fontsize',18,'XTick',-90:45:90,'xlim',[-90,90]);
title('');
xlabel('\theta, degrees');
ylabel('\rho');
print('figFolder\g694_ie4_hough_peaks.tiff','-dtiff');

%% print the hough lines located inside the cluster area
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
title('');
colorbar off;
print('figFolder\g694_ie4_hough_lines.tiff','-dtiff');

%% print branch points on skeleton
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
title('');
colorbar off;
print('figFolder\g694_ie4_skeleton_with_branch_points.tiff','-dtiff');

%% print branch points on tsN_grouped skeleton
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
title('');
caxis([12 24]);
colorbar off;
print('figFolder\g694_ie4_branchGrouped_with_branch_points.tiff','-dtiff');

%% print fragments
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
caxis([12 24]);
title('');
colorbar off;
print('figFolder\g694_ie4_fragments.tiff','-dtiff');

%% print one specific fragment, the houghline, and the gb intersection
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
caxis([0 1]);
title('');
colorbar off;
print('figFolder\g694_ie4_specific_fragment_at_gb.tiff','-dtiff');

%% print all the detected gb intersection
set(gca,'fontsize',18,'XTickLabel','','YTickLabel','');
caxis([12 24]);
title('');
colorbar off;
print('figFolder\g694_ie4_specific_all_gb_intersection.tiff','-dtiff');


