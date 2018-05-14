

% Give a specific example of g_1144. Show how dissimilarity(psi), SF(m)
% changed. Show the hcp cell.
%
% chenzhe, 2018-05-11.

close all;

% euler angle of grain-1144 in WE43_T6_C1, already corrected, so setting is not applied as '2' 
euler_d = [44.861, 164.468, 32.09];

% hcp_cell, label twin system: 1-in-[1-6], or 19-in-[19-24], or 25-in-[25-30]  
hcp_cell('euler',euler_d,'ss',25,'stress',[-1 0 0;0 0 0; 0 0 0]);
print('ts#1.tif','-dtiff');
% calculate euler after twin, and plot hcp cell again
hcp_cell('euler',euler_by_twin(euler_d,1,'Mg',2),'ss',25,'stress',[-1 0 0;0 0 0; 0 0 0],'plotTrace',0,'plotPlane',0,'plotBurgers',0);
print('ts#1_after.tif','-dtiff')

% hcp_cell, label twin system: 4-in-[1-6], or 22-in-[19-24], or 28-in-[25-30]  
hcp_cell('euler',euler_d,'ss',28,'stress',[-1 0 0;0 0 0; 0 0 0]);
print('ts#4.tif','-dtiff')
% calculate euler after twin, and plot hcp cell again
hcp_cell('euler',euler_by_twin(euler_d,4,'Mg',2),'ss',25,'stress',[-1 0 0;0 0 0; 0 0 0],'plotTrace',0,'plotPlane',0,'plotBurgers',0);
print('ts#4_after.tif','-dtiff')


%% This is a check of the identified twin system
close all;
d = [];
for iE = 2:5
    fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    load(['D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\',fName_c2t_result],'stru');
    iS = find(arrayfun(@(x) x.gID==1144, stru));
    disp(['iE = ',num2str(iE),' ------ row-0: cCentroid(1,:); row-1: twin strain; ------']);
    disp('------ row-2: pdist2 between cCentroid(1,:) and all twin strain; row-3: identified TS. ------');
    disp(stru(iS).cCen(1,:));
    disp(stru(iS).tStrain);
    disp(pdist2(stru(iS).cCen(1,:),stru(iS).tStrain));
    disp(stru(iS).trueTwin);
    d(iE,:) = pdist2(stru(iS).cCen(1,:),stru(iS).tStrain);
end

% plot the [distance] between cCentroid and twinStrain at different [iE]

d = (d(2:5,:))';
iEs = [-0.004, -0.012, -0.023, -0.039];
figure; hold on;

plot(iEs,d(1,:),'-or','linewidth',1.5,'MarkerFaceColor','r');
plot(iEs,d(2,:),'-ob','linewidth',1.5,'MarkerFaceColor','b');
plot(iEs,d(3,:),'-ok','linewidth',1.5,'MarkerFaceColor','k');
plot(iEs,d(4,:),'--sr','linewidth',1.5,'MarkerFaceColor','r');
plot(iEs,d(5,:),'--sb','linewidth',1.5,'MarkerFaceColor','b');
plot(iEs,d(6,:),'--sk','linewidth',1.5,'MarkerFaceColor','k');


xlabel('\epsilon\fontsize{14}^G');
ylabel('\psi\fontsize{14}^D');
set(gca,'xdir','reverse','ylim',[0 0.1],'xlim',[-0.053 0],'fontsize',18);
legend({'\fontsize{14}# 1','\fontsize{14}# 2','\fontsize{14}# 3','\fontsize{14}# 4','\fontsize{14}# 5','\fontsize{14}# 6'},'location','east');

print('phi vs eG for 6 TSS.tiff','-dtiff')
%% [Option B] This is a check of the identified twin system, wider plot
close all;
d = [];
for iE = 2:5
    fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    load(['D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\',fName_c2t_result],'stru');
    iS = find(arrayfun(@(x) x.gID==1144, stru));
    disp(['iE = ',num2str(iE),' ------ row-0: cCentroid(1,:); row-1: twin strain; ------']);
    disp('------ row-2: pdist2 between cCentroid(1,:) and all twin strain; row-3: identified TS. ------');
    disp(stru(iS).cCen(1,:));
    disp(stru(iS).tStrain);
    disp(pdist2(stru(iS).cCen(1,:),stru(iS).tStrain));
    disp(stru(iS).trueTwin);
    d(iE,:) = pdist2(stru(iS).cCen(1,:),stru(iS).tStrain);
end

% plot the [distance] between cCentroid and twinStrain at different [iE]

d = (d(2:5,:))';
iEs = [-0.004, -0.012, -0.023, -0.039];
figure('Position',[680,558,800,420]); hold on;

plot(iEs,d(1,:),'-or','linewidth',1.5,'MarkerFaceColor','r');
plot(iEs,d(2,:),'-ob','linewidth',1.5,'MarkerFaceColor','b');
plot(iEs,d(3,:),'-ok','linewidth',1.5,'MarkerFaceColor','k');
plot(iEs,d(4,:),'--sr','linewidth',1.5,'MarkerFaceColor','r');
plot(iEs,d(5,:),'--sb','linewidth',1.5,'MarkerFaceColor','b');
plot(iEs,d(6,:),'--sk','linewidth',1.5,'MarkerFaceColor','k');


xlabel('\epsilon\fontsize{14}^G');
ylabel('\psi\fontsize{14}^D');
set(gca,'xdir','reverse','ylim',[0 0.1],'xlim',[-0.07 0],'fontsize',18);
legend({'\fontsize{14}TS 1: (10-12)[-1011]',...
    '\fontsize{14}TS 2: (01-12)[0-111]',...
    '\fontsize{14}TS 3: (-1102)[1-101]',...
    '\fontsize{14}TS 4: (-1012)[10-11]',...
    '\fontsize{14}TS 5: (0-112)[01-11]',...
    '\fontsize{14}TS 6: (1-102)[-1101]'},'location','east');

print('phi vs eG for 6 TSS _a.tiff','-dtiff')

