
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
save('temp_results\saved_data.mat','distMap','-append','-v7.3')
%%
D = matfile('temp_results\saved_data.mat')




