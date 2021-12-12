% chenzhe 2018-08-06
% suggestted by Prof Allison
% Plot the evolution of twin area fraction vs global strain level.
% Divide the manually corrected twin map into 3x4=12 regions for errorbar
%
% Based on relabeld results, remake figure.
% chenzhe, 2019-07-27


clear twinPct;

load('D:\p\m\DIC_Analysis\20190222_1217_relabeled_result.mat','trueTwinMapCell');
load('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','strainPauses');
load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat','ID');

for iE = 2:5
    trueTwinMap = trueTwinMapCell{iE};
    if iE==2
        [nR,nC] = size(trueTwinMap);
    end
    nr = 3;
    nc = 4;
    for ir=1:nr
       for ic = 1:nc
           subTrueTwinMap = trueTwinMap([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           subIDMap = ID([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           twinPct((ir-1)*nc+ic,iE) = sum(subTrueTwinMap(:)>0)/sum(subIDMap(:)>0);
       end
    end
    
end

%%
close all;
tAvg = mean(twinPct);
tStd = std(twinPct);
eg(2:5) = strainPauses(2:5);    % eg(2:5) = [-0.004, -0.012, -0.023, -0.039]; % this is from DIC strain

errorbar(eg(2:5), 100*tAvg(2:5), 100*tStd(2:5), '-r.','linewidth',2,'markersize',18);
for ii=2:5
    text(strainPauses(ii)+0.003,0.7+tAvg(ii)*100,[num2str(100*tAvg(ii),2),'\pm',num2str(100*tStd(ii),2),'%'],'fontsize',14);
end
set(gca,'xdir','reverse','linewidth',1.5);
set(gca,'xlim',[-0.045, 0],'ylim',[0 11],'fontsize',18,'fontweight','normal');
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('Twin Area Percent, %');


%%
print('Twin Area vs Global Strain.tif','-dtiff');



%% This is a version that gives a little bit different result, as it only calculates one average of the whole map

%
load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat','ID')
load('20190222_1217_relabeled_result.mat','trueTwinMapCell')

%%

for iE=2:5
    aF(iE) = sum(trueTwinMapCell{iE}(:)>0) / sum(ID(:)>0)
end

%%
close all
figure; 
plot(strainPauses(2:5),aF(2:5)*100,'-o');
set(gca,'xdir','reverse','fontsize',16)
xlabel('strain')
ylabel('twin area percentage')
for ii=2:5
    text(strainPauses(ii)+0.003,0.5+aF(ii)*100,num2str(100*aF(ii),3));
end



