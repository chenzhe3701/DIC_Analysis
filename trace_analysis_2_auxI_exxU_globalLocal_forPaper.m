% script to make figure for paper: plot exx, and u, for strain level-5
% stitched strain and displacement field, and blow up region

strain_file = 'D:\WE43_T6_C1_insitu_compression\stitched_DIC\_5.mat'
translation_file = 'D:\WE43_T6_C1_insitu_compression\stitched_img\translations_searched_vertical_stop_0.mat'
ws_file = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\WE43_T6_C1_EbsdToSemForTraceAnalysis.mat'

addChenFunction;
load(strain_file,'exx','u');
load(translation_file);
load(ws_file,'boundaryTF','boundaryTFB','X','Y');

%% plot 
% plot exx map
[f1,a1,c1]=myplotm(exx,'x',X,'y',Y,'tf',boundaryTFB,'r',1);
clim = caxis

c = colorbar;
caxis([-0.2, 0.02])
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x, mm/mm');


% plot u map
[f2,a2,c2]=myplotm(u,'x',X,'y',Y,'tf',boundaryTFB,'r',1);
clim = caxis
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels');
caxis([-2820,310]);

%% enlarged

rTarget = 1;
cTarget = 1;
B = 1;
xmid = transX(rTarget+B,cTarget+B)+350  % overlap is about 700 pixels
ymid = transY(rTarget+B,cTarget+B)+350
indR = round((ymid-2048)/5)
indC = round((xmid-2048)/5)
xLocal = X(indR:(indR+4096/5), indC:(indC+4096/5));
yLocal = Y(indR:(indR+4096/5), indC:(indC+4096/5));
exxLocal = exx(indR:(indR+4096/5), indC:(indC+4096/5));
uLocal = u(indR:(indR+4096/5), indC:(indC+4096/5));
boundaryLocal = boundaryTFB(indR:(indR+4096/5), indC:(indC+4096/5));

% plot local exx
[f3,a3,c3]=myplotm(exxLocal,'x',xLocal,'y',yLocal,'tf',boundaryLocal,'r',1);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x, mm/mm');


% plot local u
[f4,a4,c4]=myplotm(uLocal,'x',xLocal,'y',yLocal,'tf',boundaryLocal,'r',1);
clim = caxis

c = colorbar
caxis([-185, 40]);
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels');

%% no grain boundary overlay plot 
% plot exx map
[f1,a1,c1]=myplotm(exx,'x',X,'y',Y,'r',1);
clim = caxis

c = colorbar;
caxis([-0.2, 0.02])
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x, mm/mm');


% plot u map
[f2,a2,c2]=myplotm(u,'x',X,'y',Y,'r',1);
clim = caxis
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels');
caxis([-2820,310]);

% enlarged

rTarget = 1;
cTarget = 1;
B = 1;
xmid = transX(rTarget+B,cTarget+B)+350  % overlap is about 700 pixels
ymid = transY(rTarget+B,cTarget+B)+350
indR = round((ymid-2048)/5)
indC = round((xmid-2048)/5)
xLocal = X(indR:(indR+4096/5), indC:(indC+4096/5));
yLocal = Y(indR:(indR+4096/5), indC:(indC+4096/5));
exxLocal = exx(indR:(indR+4096/5), indC:(indC+4096/5));
uLocal = u(indR:(indR+4096/5), indC:(indC+4096/5));
boundaryLocal = boundaryTFB(indR:(indR+4096/5), indC:(indC+4096/5));

% plot local exx
[f3,a3,c3]=myplotm(exxLocal,'x',xLocal,'y',yLocal,'r',1);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x, mm/mm');


% plot local u
[f4,a4,c4]=myplotm(uLocal,'x',xLocal,'y',yLocal,'r',1);
clim = caxis

c = colorbar
caxis([-185, 40]);
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels');






