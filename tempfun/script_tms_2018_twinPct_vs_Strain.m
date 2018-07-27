% plot twin size vs. compressive strain.
%
% code original for TMS 2018

close all;
figure;
plot(slevel,twinSizePct*100,'-o')
set(gca,'xlim',[0 0.04],'fontsize',18)
xlabel('Compressive Strain, mm/mm');
ylabel('Twin Area Percentage, %')