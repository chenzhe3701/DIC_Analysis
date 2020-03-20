close all;
C = colororder;

figure; % disableDefaultInteractivity(gca); 
ax1 = axes('yaxislocation','left'); hold on; 
histogram(t1.dStrain, [0:0.005:1.2],'facecolor',C(1,:));



set(ax1,'xlim',[0,1.2],'fontsize',18,'color','none');
xlabel('Diff Strain');
ylabel('Counts');

ax2 = axes('yaxislocation','right'); hold on; 
histogram(nan, [0:0.005:1.2],'facecolor',C(1,:));
histogram(t2.dStrain, [0:0.005:1.2],'facecolor',C(2,:));

set(ax2,'xlim',[0,1.2],'fontsize',18,'yaxislocation','right','color','none');
ylabel('Counts','color',C(2,:));


pos1 = get(ax1,'position');
pos2 = get(ax2,'position');
pos = [max(pos1(1),pos2(1)), max(pos1(2),pos2(2)), min(pos1(3),pos2(3)), min(pos1(4),pos2(4))];
set(ax1,'position',pos);
set(ax2,'position',pos);


axes(ax1);
breakxaxis([0.22 1.18]);
axes(ax2);
breakxaxis([0.22 1.18]);
legend('Twinned','Not twinned','location','north');hold on;
legend('Twinned','Not twinned','location','north');

