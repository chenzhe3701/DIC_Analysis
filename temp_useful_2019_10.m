
% If want to twin, (1) cannot have high basal_SF (2) cannot have low twin_SF    

close all;
iE = 5;

%% To study the effect of grain microstructure parameter, such as basal/twin_SF, grain size, etc.
%% (1) Grain twinned (red) and not twinned (blue). X: max basal SF of grain. Y: counts.
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
xlabel('Basal Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent grains twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Grains twinned & not twinned, iE = ',num2str(iE)],'fontweight','normal');

%% (2) Variants twinned & not twinned. X: basal SF. Y: coutns.
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
xlabel('Basal Schmid Factor');
ylabel('Counts');

yyaxis right;
set(gca,'ycolor','k','XTick',edges(1:end-1)+d_int,'xTickLabels',xstr,'xTickLabelRotation',45);
plot(xpos, N_t./(N_t+N_nt),'-ko','linewidth',1.5);
ylabel('Percent twinned');
legend('Not twinned','Twinned','Pct grains twinned','location','north');
set(gca,'fontsize',12,'ylim',[0 1]);
title(['Variants twinned & not twinned, iE = ',num2str(iE)],'fontweight','normal');

%% (3) Grains twinned & non-twinned. X: max twin SF. Y: counts
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

%% (4) twin SF, T2, for all variants
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






%% To study the effect of m', m'-rank, etc.
iE = 5;
load(['D:\p\m\DIC_Analysis\temp_results\twin_gb_summary_',num2str(iE),'.mat'], 'T', 'T2', 'struCell');
TT = [T;T2];

%% (I) Effect of m'-factor.  Compared <active variants> vs. <all possible variants>  

% The following generates and compares the probability mass distribution, This just means that the active variants has more with high m', and less with small m'.  

% If we want to demonstrate m' has a decisive role, it is reasonable to exclude the low SF ones, at least on the stimulation side.
% In the current study, we could first need to look at the basal

close all;

edges = [0:0.05:0.95, 1+1000*eps];
xpos = 0.025:0.05:0.975;

% (1) If look at all the active variants.  
ind = T.incoming==1;    % ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin)&(T.iiE_each_twin<=iE);    
t = T(ind,:);
% figure; histogram(t.mPrime, 0:0.05:1);
% xlabel('m'' wrt active ss/ts in neighbor');
% ylabel('Counts');
% title('for active variants at the initiating boundary','fontweight','normal');
N = histcounts(t.mPrime, edges);


% When the active mode in neighbor is basal
ind1 = ind&(ismember(T.ssn_nb,[1,2,3]));
t = T(ind1,:);
N1 = histcounts(t.mPrime, edges);
% When the active mode in neighbor is twin
ind2 = ind&(ismember(T.ssn_nb,[19:24]));
t = T(ind2,:);
N2 = histcounts(t.mPrime, edges);
figure; 
bar(xpos, [N1(:),N2(:)], 1, 'stacked');
legend({'Basal slip in neighbor','Twin in neighbor'},'location','north');
set(gca,'fontsize',16);
xlabel('m'' factor');
ylabel('Counts');
title('All active variants','fontweight','normal');

pmf = (N1+N2)./sum([N1(:);N2(:)]);  % probability mass function
pmf_1 = N1./sum(N1(:));  % probability mass function, for active variants interacting with neighbor's basal   
pmf_2 = N2./sum(N2(:));  % probability mass function, for active variants interacting with neighbor's twin  
% yyaxis right;
% plot(xpos,pmf,'-ok')

% The follwing are the baseline/background
% (2) For: all the possible variants, in all the twinned grains, look at their m' wrt all neighbor grains' MOST POSSIBLE basal         
N3 = histcounts(TT.mPrime_wrtB, edges);
figure; 
histogram(TT.mPrime_wrtB, 0:0.05:1);
xlabel('m'' factor');
ylabel('Counts');
title('m'' between all possible variants and most possible basal slip in neighbor','fontweight','normal');
pmf_b = N3./sum(N3);

% (3) For: all the possible variants, in all the twinned grains, look at their m' wrt all neighbor grains' MOST POSSIBLE twin   
N4 = histcounts(TT.mPrime_wrtT, edges);
figure; 
histogram(TT.mPrime_wrtT, 0:0.05:1);
xlabel('m'' wrt basal in neighbor');
ylabel('Counts');
title('m'' between all possible variants and most possible twin in neighbor','fontweight','normal');
pmf_t = N4./sum(N4);


% (4) Relative strength, all active wrt both basal and twin, vs. all possible wrt basal slip in neighbor     
figure; hold on;
plot(xpos, pmf, '-ko','linewidth',1.5);
plot(xpos, pmf_b, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf./pmf_b, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
xlabel('m'' factor');
ylabel('Multiples of random');
set(gca,'fontsize',16);
legend({'Active variants, ','All possible variants, basal in neighbor','Ratio'},'location','north');
title('relative strength, all vs. basal', 'fontweight','normal');

% (5) relative strength, basal slip in neighbor
figure; hold on;
plot(xpos, pmf_1, '-ko','linewidth',1.5);
plot(xpos, pmf_b, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf_1./pmf_b, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
xlabel('m'' factor');
ylabel('Multiples of random');
set(gca,'fontsize',16);
legend({'Active variants','All possible variants','Ratio'},'location','north');
title('Basal in neighbor', 'fontweight','normal');

% (6) relative strength, twin in neighbor
figure; hold on;
plot(xpos, pmf_2, '-ko','linewidth',1.5);
plot(xpos, pmf_t, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf_2./pmf_t, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
ylabel('Multiples of random');
xlabel('m'' factor');
set(gca,'fontsize',16);
legend({'Active variants','All possible variants','Ratio'},'location','north');
title('Twin in neighbor', 'fontweight','normal');

%% (II) Effect of m'-rank
close all;

edges = [0.5:1:6.5];
xpos = 1:6;

% (1) If look at all the active variants.  
ind = T.incoming==1;
t = T(ind,:);
% figure; histogram(t.rank_mPrime);
% xlabel('m''-rank wrt active ss/ts in neighbor');
% ylabel('Counts');
% title('for active variants at the initiating boundary','fontweight','normal');
N = histcounts(t.rank_mPrime, edges);

% When the active mode in neighbor is basal
ind1 = ind&(ismember(T.ssn_nb,[1,2,3]));
t = T(ind1,:);
N1 = histcounts(t.rank_mPrime, edges);
% When the active mode in neighbor is twin
ind2 = ind&(ismember(T.ssn_nb,[19:24]));
t = T(ind2,:);
N2 = histcounts(t.rank_mPrime, edges);
figure; 
bar(xpos, [N1(:),N2(:)], 1, 'stacked');
legend({'Basal slip in neighbor','Twin in neighbor'},'location','northeast');
set(gca,'fontsize',16);
xlabel('m'' rank');
ylabel('Counts');
title('All active variants','fontweight','normal');

pmf = (N1+N2)./sum([N1(:);N2(:)]);  % probability mass function
pmf_1 = N1./sum(N1(:));  % probability mass function, for active variants interacting with neighbor's basal   
pmf_2 = N2./sum(N2(:));  % probability mass function, for active variants interacting with neighbor's twin  
% yyaxis right;
% plot(xpos,pmf,'-ok')

% The follwing are the baseline/background
% (2) For: all the possible variants, in all the twinned grains, look at their m' wrt all neighbor grains' MOST POSSIBLE basal         
N3 = histcounts(TT.rank_mPrime_wrtB, edges);
figure; 
histogram(TT.rank_mPrime_wrtB, 0.5:1:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('m'' between all possible variants and most possible basal slip in neighbor','fontweight','normal');
pmf_b = N3./sum(N3);

% (3) For: all the possible variants, in all the twinned grains, look at their m' wrt all neighbor grains' MOST POSSIBLE twin   
N4 = histcounts(TT.rank_mPrime_wrtT, edges);
figure; 
histogram(TT.rank_mPrime_wrtT, 0.5:1:6.5);
xlabel('m'' rank wrt basal in neighbor');
ylabel('Counts');
title('m'' rank between all possible variants and most possible twin in neighbor','fontweight','normal');
pmf_t = N4./sum(N4);


% (4) Relative strength, all active wrt both basal and twin, vs. all possible wrt basal slip in neighbor     
figure; hold on;
plot(xpos, pmf, '-ko','linewidth',1.5);
plot(xpos, pmf_b, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf./pmf_b, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
xlabel('m'' rank');
ylabel('Multiples of random');
set(gca,'fontsize',16);
legend({'Active variants, ','All possible variants, basal in neighbor','Ratio'},'location','north');
title('relative strength, all vs. basal', 'fontweight','normal');

% (5) relative strength, basal slip in neighbor
figure; hold on;
plot(xpos, pmf_1, '-ko','linewidth',1.5);
plot(xpos, pmf_b, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf_1./pmf_b, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
xlabel('m'' rank');
ylabel('Multiples of random');
set(gca,'fontsize',16);
legend({'Active variants','All possible variants','Ratio'},'location','north');
title('Basal in neighbor', 'fontweight','normal');

% (6) relative strength, twin in neighbor
figure; hold on;
plot(xpos, pmf_2, '-ko','linewidth',1.5);
plot(xpos, pmf_t, '--kd','linewidth',1.5);
ylabel('Probability distribution');
yyaxis right;
plot(xpos, pmf_2./pmf_t, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
ylabel('Multiples of random');
xlabel('m'' factor');
set(gca,'fontsize',16);
legend({'Active variants','All possible variants','Ratio'},'location','north');
title('Twin in neighbor', 'fontweight','normal');



%% (1) all boundaries, distribution of m'-rank, is uniform, as shown before.  

%% If look at the active variants
% (2) m'-rank, boundaries with intersecting twins   
ind = (T.incoming==1);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('only intersecting twins');

% Or, just look at variants at the boundaries, the strain level at which the variant was observed at the boundary equals the strain level at which the twin was observed in grain   
% (3) boundaries with twins that initiated at this boundary
ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('twins that initiated at this boundary');


% (4) Or, just look at boundaries with twins that initiated at this boundary, and at this iE (which means newly activated at this iE)
ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin)&(T.iiE_each_twin==iE);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title(['twins that initiated at this boundary, at iE=',num2str(iE)]);

% (5) boundaries with initiating twins (determined by my assumption, only one gb in a grain is selected as the 'initiating' boundary)
ind = (T.initiating==1);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('initiating twins');
ylim = get(gca,'ylim');

% (6) divide by grain size.  Does large grain show more distinct trend? -->
% Looks like no obvious grain size effect, after manual labeling.   
ind = (T.initiating==1)&(T.gDia>100)&(T.gDia_neighbor>100);
t = T(ind,:);
N1 = histcounts(t.rank_mPrime, 0.5:6.5);
ind = (T.initiating==1)&((T.gDia<100)|(T.gDia_neighbor<100));
t = T(ind,:);
N2 = histcounts(t.rank_mPrime, 0.5:6.5);
figure;
bar(1:6,[N1(:),N2(:)], 1);
xlabel('m'' rank');
ylabel('Counts');
legend({'g diameters > 100', 'g diameters < 100'});
title('initiating twins, both grain diameter > 80um');
set(gca,'ylim',ylim);

% (7) If only look at those paired with basal slip
ind = (T.incoming==1)&(T.initiating==1)&(ismember(T.ssn_nb,[1,2,3]));
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('initiating twins not at triple points, neighbor is basal');

% (8) If only look at those paired with basal slip, and basal is big
ind = (T.incoming==1)&(T.initiating==1)&(ismember(T.ssn_nb,[1,2,3]))&(T.SF_nb>0.25);
t = T(ind,:);
figure;
histogram(t.rank_mPrime, 0.5:6.5);
xlabel('m'' rank');
ylabel('Counts');
title('initiating twins not at triple points, neighbor is basal');


%%  * likely neighbor have large max_basal_SF_nb
iE = 3;
close all
ind1 = (T.incoming==1)
t = T(ind,:);
figure;
histogram(t.max_basal_SF_nb, 0:0.05:0.5);

figure;
histogram(TT.max_basal_SF_nb, 0:0.05:0.5);

%% the effect of strain distribution: seems not distinct
close all
ind1 = (T.incoming==1)
t = T(ind,:);
figure;
histogram(t.eMean_2_nb, 0.005:0.00025: 0.03);

figure;
histogram(TT.eMean_2_nb, 0.005:0.00025: 0.03);

%% the effect of exzr_ba etc.
close all
ind1 = (T.incoming==1)
t = T(ind,:);
figure;
histogram(t.exzr_ba, 0:0.05:1);

figure;
histogram(TT.exzr_ba, 0:0.05:1);

%% the effect of exz_ba etc.
close all
ind1 = (T.incoming==1)
t = T(ind,:);
figure;
histogram(t.exz_ba, 0:0.005:0.14);


figure;
histogram(TT.exz_ba, 0:0.005:0.14);






