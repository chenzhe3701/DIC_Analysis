

%% To study the effect of m', m'-rank, etc.
iE = 5;
load(['D:\p\m\DIC_Analysis\temp_results\twin_gb_summary_',num2str(iE),'.mat'], 'T', 'T2', 'struCell');
TT = [T;T2];



%% If look at the active variants
% (2) m'-rank, boundaries with intersecting twins   
ind = (T.incoming==1);

% Or, just look at variants at the boundaries, the strain level at which the variant was observed at the boundary equals the strain level at which the twin was observed in grain   
% (3) boundaries with twins that initiated at this boundary
ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin);


% (4) Or, just look at boundaries with twins that initiated at this boundary, and at this iE (which means newly activated at this iE)
ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin)&(T.iiE_each_twin==iE);

% (5) boundaries with initiating twins (determined by my assumption, only one gb in a grain is selected as the 'initiating' boundary)
ind = (T.initiating==1);

% (6) divide by grain size.  Does large grain show more distinct trend? -->
% Looks like no obvious grain size effect, after manual labeling.   
ind = (T.initiating==1)&(T.gDia>100)&(T.gDia_neighbor>100);


% (7) If only look at those paired with basal slip
ind = (T.incoming==1)&(T.initiating==1)&(ismember(T.ssn_nb,[1,2,3]));

% (8) If only look at those paired with basal slip, and basal is big
ind = (T.incoming==1)&(T.initiating==1)&(ismember(T.ssn_nb,[1,2,3]))&(T.SF_nb>0.25);

%% Effect of tensor (strain accommodation), almost no effect on variant selection
disp('(1) Effect of tensor ===============');
close all;
edges = [0:0.013/2:0.13];
% edges = [0:0.05:1];
edges(end) = edges(end)+ 1000*eps;
xpos = edges(1:end-1) + (edges(2)-edges(1))/2;

ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin)&(T.iiE_each_twin<=iE); 
t = T(ind,:);
figure; disableDefaultInteractivity(gca);
histogram(t.exz_ba);
N = histcounts(t.exz_ba, edges);
pmf_exp = N./sum(N);

ind = (T2.TSF>0.3);   % If we select possible twin systems to be with TSF > 0.3.
tt = [T;T2(ind,:)];
figure; disableDefaultInteractivity(gca);
histogram(tt.exz_ba);
N3 = histcounts(tt.exz_ba, edges);
pmf_r = N3./sum(N3);    % theoretical and restricted

figure; disableDefaultInteractivity(gca); hold on;
plot(xpos, pmf_exp, '-ko','linewidth',1.5);
plot(xpos, pmf_r, '--kd','linewidth',1.5);
ylabel('Probability mass function');
yyaxis right;
plot(xpos, pmf_exp./pmf_r, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,2]);
ylabel('Multiples of random');
xlabel('tensor component, maybe normalized');
set(gca,'fontsize',16);
legend({'Active variants','Possible variants with SF>0.3','Ratio'},'location','north');
title('Slip/twin in neighbor', 'fontweight','normal');







%% Effect of strain, can see difference in distribution, but similar to the 9 curves which should be better categorized, no need to show.  
disp('(1) Effect of rbv ===============');
close all;
edges = [0:0.01:0.08];
edges(end) = edges(end)+ 1000*eps;
xpos = edges(1:end-1) + (edges(2)-edges(1))/2;

ind = (T.iiE_each_twin_at_this_boundary == T.iiE_each_twin)&(T.iiE_each_twin<=iE); 
t = T(ind,:);
figure; disableDefaultInteractivity(gca);
histogram(t.eMean_1_nb,edges);
N = histcounts(t.eMean_1_nb, edges);
pmf_exp = N./sum(N);

ind = (T2.TSF>0.3);   % If we select possible twin systems to be with TSF > 0.3.
tt = [T;T2(ind,:)];
figure; disableDefaultInteractivity(gca);
histogram(tt.eMean_1_nb,edges);
N3 = histcounts(tt.eMean_1_nb, edges);
pmf_r = N3./sum(N3);    % theoretical and restricted

figure; disableDefaultInteractivity(gca); hold on;
plot(xpos, pmf_exp, '-ko','linewidth',1.5);
plot(xpos, pmf_r, '--kd','linewidth',1.5);
set(gca,'ylim',get(gca,'ylim').*[0,2]);
ylabel('Probability mass function');
yyaxis right;
plot(xpos, pmf_exp./pmf_r, '-ro','linewidth',1.5);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,2]);
ylabel('Multiples of random');
xlabel('tensor component, maybe normalized');
set(gca,'fontsize',16);
legend({'Active variants','Possible variants with SF>0.3','Ratio'},'location','north');
title('Slip/twin in neighbor', 'fontweight','normal');





%% the effect of exzr_ba etc.
close all
ind1 = (T.incoming==1)
t = T(ind,:);
figure;
histogram(t.exzr_ba, 0:0.05:1);

figure;
histogram(TT.exzr_ba, 0:0.05:1);

%% The correlation between twin strength and other things.
ind = T.incoming==1;
indn = ~ind;

t = T(ind,:);
figure;
plot(t.gb_length, t.tGbVolPct, '.');






