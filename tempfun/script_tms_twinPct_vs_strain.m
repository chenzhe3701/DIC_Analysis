% script on the night before 2018 TMS
% plot twin area pct vs. global strain
% 
% Some data need to be loaded first

close all;
figure;
plot(slevel,twinSizePct*100,'-o')
set(gca,'xlim',[0 0.04],'fontsize',18)
xlabel('Compressive Strain, mm/mm');
ylabel('Twin Area Percentage, %')