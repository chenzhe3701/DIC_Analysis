close all;
figure;
plot(slevel,twinSizePct*100,'-o')
set(gca,'xlim',[0 0.04],'fontsize',18)
xlabel('Compressive Strain, mm/mm');
ylabel('Twin Area Percentage, %')