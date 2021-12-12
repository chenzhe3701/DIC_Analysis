% chenzhe 2018-08-06
% suggestted by Prof Allison
% Plot the evolution of twin area fraction vs global strain level.
% Divide the manually corrected twin map into 3x4=12 regions for errorbar

clear twinPct;
for iE = 2:5
    fileName = ['D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    load(fileName,'trueTwinMap');
    if iE==2
        [nR,nC] = size(trueTwinMap);
    end
    nr = 3;
    nc = 4;
    for ir=1:nr
       for ic = 1:nc
           subTrueTwinMap = trueTwinMap([1:floor(nR/nr)]*ir, [1:floor(nC/nc)]*ic);
           twinPct((ir-1)*nc+ic,iE) = sum(subTrueTwinMap(:)>0)/numel(subTrueTwinMap);
       end
    end
    
end

%%
close all;
tAvg = mean(twinPct);
tStd = std(twinPct);
eg(2:5) = [-0.004, -0.012, -0.023, -0.039];

errorbar(eg(2:5), 100*tAvg(2:5), 100*tStd(2:5), '-r.','linewidth',2,'markersize',18);
set(gca,'xdir','reverse','linewidth',1.5);
set(gca,'xlim',[-0.042, 0],'ylim',[0 11],'fontsize',18,'fontweight','normal');
xlabel('\fontsize{24}\epsilon\fontsize{18}^G');
ylabel('Twin Area Percent, %');
print('Twin Area vs Global Strain.tif','-dtiff');