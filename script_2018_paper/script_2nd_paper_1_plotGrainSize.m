% chenzhe, 2018-04-06 add note
% After loading data, plot the grain size histogram: Count vs Grain Size

load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\WE43_T6_C1_EbsdToSemForTraceAnalysis','gDiameter');

gDiameter(1) = 0;
figure;
histogram(gDiameter);
set(gca,'fontsize',18);
xlabel('Grain Diameter, \mum');
ylabel('Count');