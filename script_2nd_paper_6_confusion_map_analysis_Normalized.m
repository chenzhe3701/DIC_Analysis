% confusion analysis for strainScore criterion vs shapeScore criterion.
%
% chenzhe, 2018-04-09 add note
% This code, based on manual labeled results, generate a 'confusion map'
% which plot clusters with different color to indicate the TP, FP, TN cases
% in classification.
%
% chenzhe, 2018-04-30
% based on confusion analysis code version 6. I think both classifiers need
% to be re-defined, so maybe just use a new code to only look at the new
% classifiers.

clear;clc;
addChenFunction;
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

%% (0) Plot the clusterNumMap/cleaned, strainScoreMap, shapeScoreMap at a selected strain level.
iE = 4;
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile)
load(strainFile,'exx','u','v');  
fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap','scoreCF','cVolGrowthRatioMap','tProbMaxMap','mDistMap','sfMap');

%% (0.0) clusterNumMap
[f,a,c] = myplotm(clusterNumMap,'x',X,'y',Y,'tf',boundaryTFB,'r',2);
colors = lines(7);
colorMap = [0 0 0; colors([3,1,5,2,4],:)];
colormap(colorMap);
caxis([-0.5 5.5]);
set(c,'Limits',[0.5 5.5],'Ticks',[1:5]);
% set(a,'fontsize',24,'XTick',[],'yTick',[]);   % old scripts
title('Cluster Label','fontweight','normal');

%% [***] maximize plot and run this:
script_make_double_axis;
title(c,'Cluster ID','fontsize',24);
print('clusterID iE=4.tif','-dtiff');
%%
% plot strain/u map, and enlarged part of a grain 1144
% get data for enlarged grain
ind_local = ismember(ID, 1144); %ismember(ID, [ID_current,ID_neighbor]);
indC_min = find(sum(ind_local, 1), 1, 'first')-50;
indC_max = find(sum(ind_local, 1), 1, 'last')+50;
indR_min = find(sum(ind_local, 2), 1, 'first')-25;
indR_max = find(sum(ind_local, 2), 1, 'last')+25;

xLocal = X(indR_min:indR_max, indC_min:indC_max);
yLocal = Y(indR_min:indR_max, indC_min:indC_max);
exxLocal = exx(indR_min:indR_max, indC_min:indC_max);
uLocal = u(indR_min:indR_max, indC_min:indC_max);
boundaryLocal = boundaryTFB(indR_min:indR_max, indC_min:indC_max);

% (0.1) plot local exx
[f,a,c]=myplotm(exxLocal,'x',xLocal,'y',yLocal,'tf',boundaryLocal,'r',1);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x','fontweight','normal');


% (0.2) plot local u
[f,a,c]=myplotm(uLocal,'x',xLocal,'y',yLocal,'tf',boundaryLocal,'r',1);
caxis_m([-1540, -1440]);

set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels','fontweight','normal');

%% (0.3) plot exx map whole, with grain boundary overlay
[f,a,c]=myplotm(exx,'x',X,'y',Y,'tf',boundaryTFB,'r',1);
c = caxis_m([-0.12, 0.02]);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('\epsilon_x_x','fontweight','normal');
rectangle('Position',[xLocal(1), yLocal(1), xLocal(end)-xLocal(1), yLocal(end)-yLocal(1)],'linewidth',4)
%% [***] maximize plot and run this:
script_make_double_axis;
title(c,'\epsilon_{xx}','fontsize',24);
print('exx iE=4.tif','-dtiff');

%% (0.4) plot u map whole, with grain boundary overlay
[f,a,c]=myplotm(u,'x',X,'y',Y,'tf',boundaryTFB,'r',1);
c = caxis_m([-1800,200]);
set(gca,'fontsize',24)
xlabel('X, pixels');
ylabel('Y, pixels');
title('u, pixels','fontweight','normal');
rectangle('Position',[xLocal(1), yLocal(1), xLocal(end)-xLocal(1), yLocal(end)-yLocal(1)],'linewidth',4)
%% [***] maximize plot and run this:
script_make_double_axis;
title(c,'u, pixels','fontsize',24);
print('u iE=4.tif','-dtiff');


%% (0.5) plot ShapeScoreMapNormalized, adjust color so that high vs low scores are distinct
shapeScoreMapNew = (cVolGrowthRatioMap+1)/2.*tProbMaxMap;
[f,a,c] = myplot(X,Y,shapeScoreMapNew,boundaryTFB,3);
% set(gca,'fontsize',24,'XTick',[],'YTick',[])
title('ShapeScore','fontweight','normal');
caxis([0 1]);
% colormap(parula(8));
%% maximize plot and run this:
script_make_double_axis;
title(c,'\eta','fontsize',24);
print('shape score Normalized iE=4.tif','-dtiff');


%% (0.6) plot strainScoreMapNormalized, adjust color so that high vs low scores are distinct
C = 7; H = 1;
strainScoreMapNew = sfMap/(1+C*H) -C/(1+C*H)*min(H,mDistMap) + (1+2*C*H)/2/(1+C*H);
[f,a,c] = myplot(X,Y,strainScoreMapNew,boundaryTFB,3);
% set(gca,'fontsize',24,'XTick',[],'YTick',[])
title('ShapeScore','fontweight','normal');
caxis([0.9 1]);
% colormap(parula(8));
%% maximize plot and run this:
script_make_double_axis;
title(c,'\eta','fontsize',24);
print('strain score Normalized iE=4.tif','-dtiff');


%% (1) Generate a confusion map for selected iE. Label TP, FP, TN with different color.
% This represents the 'forward' method: set-up a threshold by visual
% observation, and we know that qualitatively, it works a little bit.
% But to fully label all the cases, we have to guided by these criterions,
% and do a lot of manual clean-up.
% iE = 4 is a good example. (I don't think we should provide such maps at
% other iE levels

iE = 4;

colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap','scoreCF','cVolGrowthRatioMap','tProbMaxMap','mDistMap','sfMap');

%% (1.1) Use [StrainScore_Normalized (c=7, H=1)], set (phi_th=0.965)
strainScoreCF = 0.965;    % 0.22, 0.17
shapeScoreCF = 0.2;
% previous criterion, 7Dis+(0.5-SF)< scoreCF
% i.e., new criterion, 7Dis-SF<socre-0.5

TP = 0;    % hit
FP = 0;    % false alarm
FN = 0;    % missing
TN = 0;    % correct rejection

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end
C = 7;      % default = 7, but we can look at other score: 5.03; 6.36; 6.24; 3.59; for iE = 2,3,4,5
H = 1;

for iS=1:length(stru)
    for iCluster=1:length(stru(iS).cLabel)
        strainScore = stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H);
        shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
        if stru(iS).trueTwin(iCluster)>0    % if ground true
            if strainScore > strainScoreCF
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else        % ground false
            if strainScore > strainScoreCF
                FP = FP + 1;
            else
                TN = TN + 1;
            end
        end
        
    end
end

% plot 'confusion map'
confMap = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
confMap(trueTwinMap>0)=1;   % true twin

% add by 2. So if it's TP, it will be 1+2=3. If it's FP, it will be 0+2=2.
% (1) strainScoreNormalized as criterion
strainScoreMapNew = sfMap/(1+C*H) -C/(1+C*H)*min(H,mDistMap) + (1+2*C*H)/2/(1+C*H);
confMap(strainScoreMapNew>strainScoreCF) = confMap(strainScoreMapNew>strainScoreCF)+2;

% make only boundary as 0, others nan
confMap(confMap==0) = 0;
confMap(boundaryTFB==1) = nan;

[f,a,c] = myplot(X, Y, confMap, boundaryTFB, 3);
colormap(colorMap);
caxis([-1.5 3.5]);set(c,'limits',[-0.5, 3.5]);
c.Ticks = 0:3;
c.TickLabels={['TN: ',num2str(TN)], ['FN: ',num2str(FN)],['FP: ',num2str(FP)],['TP: ',num2str(TP)]};

% set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
set(a,'fontsize',18);

ttl = ['StrainScoreNormalized-th=',num2str(strainScoreCF)];
ttl_text = ['\phi_{th}= ',num2str(strainScoreCF)];
title(a,ttl,'fontweight','normal');
annotation(f,'textbox', [0.635 0.96 0.3 0.042], 'String',ttl_text, 'LineStyle','none', 'FontSize',24);
%% maximize plot and run this:
script_make_double_axis;
print([ttl,'.tif'],'-dtiff');



%% (1.2) Use [StrainScore_Normalized (c=7, H=1)], set (phi_th=0.958)
strainScoreCF = 0.958;    % 0.22, 0.17
shapeScoreCF = 0.2;
% previous criterion, 7Dis+(0.5-SF)< scoreCF
% i.e., new criterion, 7Dis-SF<socre-0.5

TP = 0;    % hit
FP = 0;    % false alarm
FN = 0;    % missing
TN = 0;    % correct rejection

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end
C = 7;      % default = 7, but we can look at other score: 5.03; 6.36; 6.24; 3.59; for iE = 2,3,4,5
H = 1;

for iS=1:length(stru)
    for iCluster=1:length(stru(iS).cLabel)
        strainScore = stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H);
        shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
        if stru(iS).trueTwin(iCluster)>0    % if ground true
            if strainScore > strainScoreCF
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else        % ground false
            if strainScore > strainScoreCF
                FP = FP + 1;
            else
                TN = TN + 1;
            end
        end
        
    end
end

% plot 'confusion map'
confMap = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
confMap(trueTwinMap>0)=1;   % true twin

% add by 2. So if it's TP, it will be 1+2=3. If it's FP, it will be 0+2=2.
% (1) strainScoreNormalized as criterion
strainScoreMapNew = sfMap/(1+C*H) -C/(1+C*H)*min(H,mDistMap) + (1+2*C*H)/2/(1+C*H);
confMap(strainScoreMapNew>strainScoreCF) = confMap(strainScoreMapNew>strainScoreCF)+2;

% make only boundary as 0, others nan
confMap(confMap==0) = 0;
confMap(boundaryTFB==1) = nan;

[f,a,c] = myplot(X, Y, confMap, boundaryTFB, 3);
colormap(colorMap);
caxis([-1.5 3.5]);set(c,'limits',[-0.5, 3.5]);
c.Ticks = 0:3;
c.TickLabels={['TN: ',num2str(TN)], ['FN: ',num2str(FN)],['FP: ',num2str(FP)],['TP: ',num2str(TP)]};

% set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
set(a,'fontsize',18);

ttl = ['StrainScoreNormalized-th=',num2str(strainScoreCF)];
ttl_text = ['\phi_{th}= ',num2str(strainScoreCF)];
title(a,ttl,'fontweight','normal');
annotation(f,'textbox', [0.635 0.96 0.3 0.042], 'String',ttl_text, 'LineStyle','none', 'FontSize',24);
%% maximize plot and run this:
script_make_double_axis;
print([ttl,'.tif'],'-dtiff');



%% (1.3) use [eta_ShapeScore_Normalized], set (eta_th=0.55)
strainScoreCF = 0.965;    % 0.22, 0.17
shapeScoreCF = 0.55;

TP = 0;    % hit
FP = 0;    % false alarm
FN = 0;    % missing
TN = 0;    % correct rejection

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end
C = 7;      % default = 7, but we can look at other score: 5.03; 6.36; 6.24; 3.59; for iE = 2,3,4,5
H = 1;

for iS=1:length(stru)
    for iCluster=1:length(stru(iS).cLabel)
        strainScore = stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H);
        shapeScore = (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster);
        if stru(iS).trueTwin(iCluster)>0    % if ground true
            if shapeScore > shapeScoreCF
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else        % ground false
            if shapeScore > shapeScoreCF
                FP = FP + 1;
            else
                TN = TN + 1;
            end
        end
    end
end

% plot 'confusion map'
confMap = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
confMap(trueTwinMap>0)=1;   % true twin

% add by 2. So if it's TP, it will be 1+2=3. If it's FP, it will be 0+2=2.
% (3) shapeScoreNormalized (New) as criterion
shapeScoreMapNew = (cVolGrowthRatioMap+1)/2.*tProbMaxMap;
confMap(shapeScoreMapNew>shapeScoreCF) = confMap(shapeScoreMapNew>shapeScoreCF)+2;

% make only boundary as 0, others nan
confMap(confMap==0) = 0;
confMap(boundaryTFB==1) = nan;

[f,a,c] = myplot(X, Y, confMap, boundaryTFB, 3);
colormap(colorMap);
caxis([-1.5 3.5]);set(c,'limits',[-0.5, 3.5]);
c.Ticks = 0:3;
c.TickLabels={['TN: ',num2str(TN)], ['FN: ',num2str(FN)],['FP: ',num2str(FP)],['TP: ',num2str(TP)]};

% set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
set(a,'fontsize',18);
ttl = ['ShapeScoreNormalized-th=',num2str(shapeScoreCF)];
ttl_text = ['\eta_{th}= ',num2str(shapeScoreCF)];
title(a,ttl,'fontweight','normal');
annotation(f,'textbox', [0.635 0.96 0.3 0.042], 'String',ttl_text, 'LineStyle','none', 'FontSize',24);

%% maximize plot and run this:
script_make_double_axis;
print([ttl,'.tif'],'-dtiff');




%% (1.4) use both score [StrainScoreNormalized (c=7 H=1)] OR [ShapeScore_New], set (phi_th=-0.22)OR(eta_th=0.55)
strainScoreCF = 0.965;    % 0.22, 0.17
shapeScoreCF = 0.55;
% previous criterion, 7Dis+(0.5-SF)< scoreCF
% i.e., new criterion, 7Dis-SF<socre-0.5

TP = 0;    % hit
FP = 0;    % false alarm
FN = 0;    % missing
TN = 0;    % correct rejection

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end
C = 7;      % default = 7, but we can look at other score: 5.03; 6.36; 6.24; 3.59; for iE = 2,3,4,5
H = 1;

for iS=1:length(stru)
    for iCluster=1:length(stru(iS).cLabel)
        strainScore = stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H);
        shapeScore = (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster);
        if stru(iS).trueTwin(iCluster)>0    % if ground true
            if (strainScore > strainScoreCF)||(shapeScore > shapeScoreCF)
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else        % ground false
            if (strainScore > strainScoreCF)||(shapeScore > shapeScoreCF)
                FP = FP + 1;
            else
                TN = TN + 1;
            end
        end
    end
end

% plot 'confusion map'
confMap = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
confMap(trueTwinMap>0)=1;   % true twin

% add by 2. So if it's TP, it will be 1+2=3. If it's FP, it will be 0+2=2.
% (4) combine both as criterion, strainScore_c=7 OR shapeScore_new
strainScoreMapNew = sfMap/(1+C*H) -C/(1+C*H)*min(H,mDistMap) + (1+2*C*H)/2/(1+C*H);
shapeScoreMapNew = (cVolGrowthRatioMap+1)/2.*tProbMaxMap;
confMap((strainScoreMapNew>strainScoreCF)|(shapeScoreMapNew>shapeScoreCF)) = confMap((strainScoreMapNew>strainScoreCF)|(shapeScoreMapNew>shapeScoreCF))+2;

% make only boundary as 0, others nan
confMap(confMap==0) = 0;
confMap(boundaryTFB==1) = nan;

[f,a,c] = myplot(X, Y, confMap, boundaryTFB, 3);
colormap(colorMap);
caxis([-1.5 3.5]);set(c,'limits',[-0.5, 3.5]);
c.Ticks = 0:3;
c.TickLabels={['TN: ',num2str(TN)], ['FN: ',num2str(FN)],['FP: ',num2str(FP)],['TP: ',num2str(TP)]};

% set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
set(a,'fontsize',18);
ttl = ['(StrainScoreNormalized-th=',num2str(strainScoreCF),') OR (ShapeScoreNormalized-th=',num2str(shapeScoreCF),')'];
ttl_text = ['(\phi_{th}= ',num2str(strainScoreCF),') OR (\eta_{th}= ',num2str(shapeScoreCF),')'];
title(a,ttl,'fontweight','normal');
annotation(f,'textbox', [0.635 0.96 0.3 0.042], 'String',ttl_text, 'LineStyle','none', 'FontSize',24);
%% maximize plot and run this:
script_make_double_axis;
print([ttl,'.tif'],'-dtiff');



%% (2) summarize the twin area/vol fraction at strain levels 2-5
for iE = 2:5
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap');
    twinSizePct(iE) = sum(trueTwinMap(:)>0)/sum(trueTwinMap(:)>=0);
end

figure;
plot(strainPauses(2:5),twinSizePct(2:5)*100,'-or','linewidth',1.5);
set(gca,'xdir','reverse','fontsize',18);
xlabel('\fontsize{24}\epsilon\fontsize{16}^G'); % xlabel('Global Uniaxial Strain, mm/mm');
ylabel('Twin Area Percent, %');
%%
print(['twin pct vs strain.tif'],'-dtiff');


%% (3) Confusion analysis using ACC metric, plot ACC curve vs criterion threshold
useACC = 1;
useF1 = 1-useACC;

%% (3.1) [ACC] vs. [StrainScoreNormalized_th] alone (C=7) , compare at [ALL] strain levels
close all;
clear scoreSorted ACC ACC_max F1 F1_max bestScore TP_best FP_best FN_best TN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);
colors([2,3],:)=colors([3,2],:);
colors(4,:) = [0 0 1];
colors(5,:) = [0 0 0];    % black

for iE=iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    
    switch iE
        case 2
            C = 5.03;
        case 3
            C = 6.36;
        case 4
            C = 6.24;
        case 5
            C = 3.59;
    end
    C = 7;
    H = 1;
    % create vectors containing necessary info
    info = [];
    for iS = 1:length(stru)
        for iCluster = 1:length(stru(iS).cLabel)
            % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
            info = [info; stru(iS).trueTwin(iCluster),...
                stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
                (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
                1];
        end
    end
    
    % sort rows, based on col-2(strainScoreNormalized(C=7,H))    
    info = sortrows(info,2,'descend');
    scoreSorted{iE} = info(:,2);
    
    % true positive, ..., etc
    TP = cumsum(info(:,1)>0);   % col-2 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives.
    P = TP(end);                % The end of the cumsum is the total positive cases
    FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.
    
    FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives.
    N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
    TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.
    
    % rate
    TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
    FNR = FN./P;
    PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
    
    % accuracy
    ACC{iE} = (TP+TN)./(TP+FP+FN+TN);
    % F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
    % sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
    F1{iE} = 2*TP./(2*TP+FP+FN);
    
    [ACC_max(iE),ind] = max(ACC{iE});
    bestScore(iE) = scoreSorted{iE}(ind);
    TP_best(iE) = TP(ind);
    FP_best(iE) = FP(ind);
    FN_best(iE) = FN(ind);
    TN_best(iE) = TN(ind);
    PPV_best(iE) = PPV(ind);
    TPR_best(iE) = TPR(ind);
    
    % plot(scoreSorted, F1,'linewidth',1.5,'color',colors(iE-1,:));
    
end

plot(scoreSorted{2}, ACC{2},'-', 'linewidth',2,'color',colors(2,:));
plot(scoreSorted{3}, ACC{3},'--','linewidth',2,'color',colors(3,:));
plot(scoreSorted{4}, ACC{4},'-.','linewidth',2,'color',colors(4,:));
plot(scoreSorted{5}, ACC{5},':', 'linewidth',2,'color',colors(5,:));

ylabel('ACC');
xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
set(gca,'fontsize',18,'xlim',[0,1]);
legend({'strain level 2: -0.6%','strain level 3: -1.2%','strain level 4: -2.1%','strain level 5: -3.7%'},'location','best');
legend({'\fontsize{24}\epsilon\fontsize{16}^G = -0.004','\fontsize{24}\epsilon\fontsize{16}^G = -0.012','\fontsize{24}\epsilon\fontsize{16}^G = -0.023','\fontsize{24}\epsilon\fontsize{16}^G = -0.039'},'location','best');

t1 = table(bestScore(iE_start:iE_stop)',...
    ACC_max(iE_start:iE_stop)',...
    TP_best(iE_start:iE_stop)',...
    FP_best(iE_start:iE_stop)',...
    FN_best(iE_start:iE_stop)',...
    TN_best(iE_start:iE_stop)',...
    PPV_best(iE_start:iE_stop)',...
    TPR_best(iE_start:iE_stop)',...
    'VariableNames',{'StrainScore_th','Accuracy','TP','FP','FN','TN','Precision','Sensitivity'})
%%
print(['Event Proportion vs StrainScoreNormalized_th c_7.tif'],'-dtiff');


%% (3.2) [ACC] vs. [ShapeScoreNormalized_th] alone, compare at [ALL] strain levels
clear scoreSorted ACC ACC_max F1 F1_max bestScore TP_best FP_best FN_best TN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);
colors([2,3],:)=colors([3,2],:);
colors(4,:) = [0 0 1];
colors(5,:) = [0 0 0];    % black

for iE=iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    
    switch iE
        case 2
            C = 5.03;
        case 3
            C = 6.36;
        case 4
            C = 6.24;
        case 5
            C = 3.59;
    end
    C=7;
    H = 1;
    % create vectors containing necessary info
    info = [];
    for iS = 1:length(stru)
        for iCluster = 1:length(stru(iS).cLabel)
            % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
            info = [info; stru(iS).trueTwin(iCluster),...
                stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
                (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
                1];
        end
    end
    
    
    % sort rows, based on col-3(shapeScore_Normalized)    
    info = sortrows(info,3,'descend');
    scoreSorted{iE} = info(:,3);
    
    % true positive, ..., etc
    TP = cumsum(info(:,1)>0);   % col-2 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives.
    P = TP(end);                % The end of the cumsum is the total positive cases
    FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.
    
    FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives.
    N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
    TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.
    
    % rate
    TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
    FNR = FN./P;
    PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
    
    % accuracy
    ACC{iE} = (TP+TN)./(TP+FP+FN+TN);
    % F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
    % sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
    F1{iE} = 2*TP./(2*TP+FP+FN);

    [ACC_max(iE),ind] = max(ACC{iE});
    bestScore(iE) = scoreSorted{iE}(ind);
    TP_best(iE) = TP(ind);
    FP_best(iE) = FP(ind);
    FN_best(iE) = FN(ind);
    TN_best(iE) = TN(ind);
    PPV_best(iE) = PPV(ind);
    TPR_best(iE) = TPR(ind);

    % plot(scoreSorted, F1,'linewidth',1.5,'color',colors(iE-1,:));
end


plot(scoreSorted{2}, ACC{2},'-', 'linewidth',2,'color',colors(2,:));
plot(scoreSorted{3}, ACC{3},'--','linewidth',2,'color',colors(3,:));
plot(scoreSorted{4}, ACC{4},'-.','linewidth',2,'color',colors(4,:));
plot(scoreSorted{5}, ACC{5},':', 'linewidth',2,'color',colors(5,:));

ylabel('ACC');
xlabel('\eta_{th}'); % xlabel('ShapeScore Threshold');
set(gca,'fontsize',18,'xlim',[0,1]);
legend({'strain level 2: -0.6%','strain level 3: -1.2%','strain level 4: -2.1%','strain level 5: -3.7%'},'location','best');
legend({'\fontsize{24}\epsilon\fontsize{16}^G = -0.004','\fontsize{24}\epsilon\fontsize{16}^G = -0.012','\fontsize{24}\epsilon\fontsize{16}^G = -0.023','\fontsize{24}\epsilon\fontsize{16}^G = -0.039'},'location','best');

t1 = table(bestScore(iE_start:iE_stop)',...
    ACC_max(iE_start:iE_stop)',...
    TP_best(iE_start:iE_stop)',...
    FP_best(iE_start:iE_stop)',...
    FN_best(iE_start:iE_stop)',...
    TN_best(iE_start:iE_stop)',...
    PPV_best(iE_start:iE_stop)',...
    TPR_best(iE_start:iE_stop)',...
    'VariableNames',{'ShapeScore_th','Accuracy','TP','FP','FN','TN','Precision','Sensitivity'})
%%
print(['Event Proportion vs ShapeScoreNormalized_th.tif'],'-dtiff');

%% (3.3) [ACC] vs. [StrainScoreNormalized_th] alone (C=Fitted) , compare at [ALL] strain levels
clear scoreSorted ACC ACC_max F1 F1_max bestScore TP_best FP_best FN_best TN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);
colors([2,3],:)=colors([3,2],:);
colors(4,:) = [0 0 1];
colors(5,:) = [0 0 0];    % black

for iE=iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    
    switch iE
        case 2
            C = 5.03;
        case 3
            C = 6.36;
        case 4
            C = 6.24;
        case 5
            C = 3.59;
    end

    H = 1;
    % create vectors containing necessary info
    info = [];
    for iS = 1:length(stru)
        for iCluster = 1:length(stru(iS).cLabel)
            % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
            info = [info; stru(iS).trueTwin(iCluster),...
                stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
                (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
                1];
        end
    end
    
    
    % sort rows, based on col-2(strainScore_Normalized,c=fitted7)    
    info = sortrows(info,2,'descend');
    scoreSorted{iE} = info(:,2);
    
    % true positive, ..., etc
    TP = cumsum(info(:,1)>0);   % col-2 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives.
    P = TP(end);                % The end of the cumsum is the total positive cases
    FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.
    
    FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives.
    N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
    TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.
    
    % rate
    TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
    FNR = FN./P;
    PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
    
    % accuracy
    ACC{iE} = (TP+TN)./(TP+FP+FN+TN);
    % F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
    % sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
    F1{iE} = 2*TP./(2*TP+FP+FN);
    
    [ACC_max(iE),ind] = max(ACC{iE});
    bestScore(iE) = scoreSorted{iE}(ind);
    TP_best(iE) = TP(ind);
    FP_best(iE) = FP(ind);
    FN_best(iE) = FN(ind);
    TN_best(iE) = TN(ind);
    PPV_best(iE) = PPV(ind);
    TPR_best(iE) = TPR(ind);
    
    % plot(scoreSorted, F1,'linewidth',1.5,'color',colors(iE-1,:));
    
end

plot(scoreSorted{2}, ACC{2},'-', 'linewidth',2,'color',colors(2,:));
plot(scoreSorted{3}, ACC{3},'--','linewidth',2,'color',colors(3,:));
plot(scoreSorted{4}, ACC{4},'-.','linewidth',2,'color',colors(4,:));
plot(scoreSorted{5}, ACC{5},':', 'linewidth',2,'color',colors(5,:));

ylabel('ACC');
xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
set(gca,'fontsize',18,'xlim',[0,1]);
legend({'strain level 2: -0.6%','strain level 3: -1.2%','strain level 4: -2.1%','strain level 5: -3.7%'},'location','best');
legend({'\fontsize{24}\epsilon\fontsize{16}^G = -0.004','\fontsize{24}\epsilon\fontsize{16}^G = -0.012','\fontsize{24}\epsilon\fontsize{16}^G = -0.023','\fontsize{24}\epsilon\fontsize{16}^G = -0.039'},'location','best');

t1 = table(bestScore(iE_start:iE_stop)',...
    ACC_max(iE_start:iE_stop)',...
    TP_best(iE_start:iE_stop)',...
    FP_best(iE_start:iE_stop)',...
    FN_best(iE_start:iE_stop)',...
    TN_best(iE_start:iE_stop)',...
    PPV_best(iE_start:iE_stop)',...
    TPR_best(iE_start:iE_stop)',...
    'VariableNames',{'StrainScore_th','Accuracy','TP','FP','FN','TN','Precision','Sensitivity'})
%%
print(['Event Proportion vs StrainScoreNormalized_th C_fitted.tif'],'-dtiff');


%% (4) For selected iE, Look at metrics change with threshoulds
close all;
%% (4.1 a)  Just look at [iE=4], how [Counts] (TP,FP,FN,TN) changes with [StrainScoreNormalized_th] (c=7)   
clear scoreSorted ACC ACC_max F1 F1_max bestScore TP_best FP_best FN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);

iE = 4;
fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru');

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end
C = 7;
H = 1;
% create vectors containing necessary info
info = [];
for iS = 1:length(stru)
    for iCluster = 1:length(stru(iS).cLabel)
        % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
        info = [info; stru(iS).trueTwin(iCluster),...
            stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
            (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
            1];
    end
end

% sort rows, based on col-2(StrainScoreNormalized, c=7)
info = sortrows(info,2,'descend');
scoreSorted = info(:,2);

% true positive, ..., etc
TP = cumsum(info(:,1)>0);   % col-3 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives. 
P = TP(end);                % The end of the cumsum is the total positive cases 
FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.

FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives. 
N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.  

% rate
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
FNR = FN./P;
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).

% accuracy
ACC = (TP+TN)./(TP+FP+FN+TN);
% F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
% sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
F1 = 2*TP./(2*TP+FP+FN);

plot(scoreSorted, TP,'linewidth',2,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',2,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',2,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
legend({'TP','FP','FN','TN'},'location','best');

xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
ylabel('Count');
set(gca,'fontsize',18,'xlim',[0 1])%,'YScale','log','YTick',[0 10 100 1000 10000]);

rectangle('Position',[0.9 0 0.1 600],'edgecolor','k','linewidth',2)
%%
print(['iE=',num2str(iE),' Event Count vs StrainScoreNormalized_th.tif'],'-dtiff');


%% (4.1 b) zoom in
figure; hold on;
plot(scoreSorted, TP,'linewidth',4,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',4,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',4,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
set(gca,'xlim',[0.9, 1],'ylim',[0 600],'fontsize',18)
%%
print(['iE=',num2str(iE),' Event Count vs StrainScoreNormalized_th zoom.tif'],'-dtiff');


%% (4.1 c) at [iE=4] how [StrainScoreNormalized_th] affect [ACC,precision,sentivity], precision = PPV = TP/TP+FP, and sensitivity TPR=TP/TP+FN   
figure; hold on;
plot(scoreSorted, ACC,'linewidth',3,'color',[0 0 0]);
plot(scoreSorted, PPV,'linewidth',1.5,'color',[1 0 0]);
plot(scoreSorted, TPR,'linewidth',1.5,'color',[0 0 1]);

xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
ylabel('Rate');
set(gca,'fontsize',18,'xlim',[0 1],'ylim',[0,1.1]);

legend({'Accuracy (ACC)','Precision (PPV)','Sensitivity (TPR)'},'location','best');
%%
print(['iE=',num2str(iE),' Event Proportion vs StrainScoreNormalized_th.tif'],'-dtiff');


%% (4.2 a)  Just look at [iE=4], how [Counts] (TP,FP,FN,TN) changes with [StrainScoreNormalized_th] (c=fitted) 
clear scoreSorted ACC ACC_max F1 F1_max bestScore TP_best FP_best FN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);

iE = 4;
fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru');

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end

H = 1;
% create vectors containing necessary info
info = [];
for iS = 1:length(stru)
    for iCluster = 1:length(stru(iS).cLabel)
        % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
        info = [info; stru(iS).trueTwin(iCluster),...
            stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
            (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
            1];
    end
end

% sort rows, based on col-2(StrainScoreNormalized, c=fitted)
info = sortrows(info,2,'descend');
scoreSorted = info(:,2);

% true positive, ..., etc
TP = cumsum(info(:,1)>0);   % col-3 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives. 
P = TP(end);                % The end of the cumsum is the total positive cases 
FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.

FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives. 
N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.  

% rate
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
FNR = FN./P;
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).

% accuracy
ACC = (TP+TN)./(TP+FP+FN+TN);
% F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
% sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
F1 = 2*TP./(2*TP+FP+FN);

plot(scoreSorted, TP,'linewidth',2,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',2,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',2,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
legend({'TP','FP','FN','TN'},'location','best');

xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
ylabel('Count');
set(gca,'fontsize',18,'xlim',[0 1])%,'YScale','log','YTick',[0 10 100 1000 10000]);

rectangle('Position',[0.9 0 0.1 600],'edgecolor','k','linewidth',2)
%%
print(['iE=',num2str(iE),' Event Count vs StrainScoreNormalized_th C_fitted.tif'],'-dtiff');


%% (4.2 b) zoom in
figure; hold on;
plot(scoreSorted, TP,'linewidth',4,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',4,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',4,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
set(gca,'xlim',[0.9 1],'ylim',[0 600],'fontsize',18)
%%
print(['iE=',num2str(iE),' Event Count vs StrainScoreNormalized_th C_fitted zoom.tif'],'-dtiff');


%% (4.2 c) how [StrainScoreNormalized_th] affect [ACC,precision,sensitivity], precision = PPV = TP/TP+FP, and sensitivity TPR=TP/TP+FN   
figure; hold on;
plot(scoreSorted, ACC,'linewidth',3,'color',[0 0 0]);
plot(scoreSorted, PPV,'linewidth',1.5,'color',[1 0 0]);
plot(scoreSorted, TPR,'linewidth',1.5,'color',[0 0 1]);

xlabel('\phi_{th}'); % xlabel('StrainScore Threshold');
ylabel('Rate');
set(gca,'fontsize',18,'xlim',[0 1],'ylim',[0,1.1]);

legend({'Accuracy (ACC)','Precision (PPV)','Sensitivity (TPR)'},'location','best');
%%
print(['iE=',num2str(iE),' Event Proportion vs StrainScoreNormalized_th C_fitted.tif'],'-dtiff');


%% (4.3 a)  Just look at [iE=4], how [Counts] (TP,FP,FN,TN) changes with [ShapeScoreNormalized_th]   
clear scoreSorted F1 F1_max bestScore TP_best FP_best FN_best PPV_best TPR_best;
figure; hold on;
colors = lines(7);

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru');

switch iE
    case 2
        C = 5.03;
    case 3
        C = 6.36;
    case 4
        C = 6.24;
    case 5
        C = 3.59;
end

H = 1;
% create vectors containing necessary info
info = [];
for iS = 1:length(stru)
    for iCluster = 1:length(stru(iS).cLabel)
        % [trueTwin, strainScoreNormalized(C,H), shapeScoreNormalized,  dummy=1, ...]
        info = [info; stru(iS).trueTwin(iCluster),...
            stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
            (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
            1];
    end
end

% sort rows, based on col-3(shapeScoreNormalized)    
info = sortrows(info,3,'descend');
scoreSorted = info(:,3);

% true positive, ..., etc
TP = cumsum(info(:,1)>0);   % col-3 as threshold. Cumsum of entries to this point are the identified positives. Among which values in col-1 and >0 are the true positives. 
P = TP(end);                % The end of the cumsum is the total positive cases 
FN = P-TP;                  % Cumsum of entries below this point are the identified negatives. Among which the remaining (P-TP) are identified negative but actually positive, i.e., false negative.

FP = cumsum(info(:,1)==0);  % Similar to 'TP', To this point, in col-1, values <= 0 are the false positives. 
N = FP(end);                % When threshould at the end, all are identified positive. The false part is the true negatives.
TN = N-FP;                  % Above the threshold, the false positives are actually negatives.  The remaining (N-FP) below the threshold are the true negatives.  

% rate
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
FNR = FN./P;
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).

% accuracy
ACC = (TP+TN)./(TP+FP+FN+TN);
% F1 score. F1 is the harmonic average of the precision=TP/(TP+FP), and
% sensitivity=recall=TPR=TP/(TP+FN), and F1 = 2/(1/precision+1/recall)
F1 = 2*TP./(2*TP+FP+FN);

plot(scoreSorted, TP,'linewidth',2,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',2,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',2,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
legend({'TP','FP','FN','TN'},'location','best');


xlabel('\eta_{th}'); % xlabel('ShapeScore Threshold');
ylabel('Count');
set(gca,'fontsize',18,'xlim',[0 1]);%,'YScale','log','YTick',[0 10 100 1000 10000]);

rectangle('Position',[0.45 0 0.3 600],'edgecolor','k','linewidth',2)
%%
print(['iE=',num2str(iE),' Event Count vs ShapeScoreNormalized_th.tif'],'-dtiff');


%% (4.3 b) zoom in
figure; hold on;
plot(scoreSorted, TP,'linewidth',4,'color',[1 0 0]);
plot(scoreSorted, FP,'linewidth',4,'color',colors(5,:));
plot(scoreSorted, FN,'linewidth',4,'color',[0 0 1]);
plot(scoreSorted, TN,'--k','linewidth',1.5);
set(gca,'fontsize',18,'xlim',[0.45, 0.75],'ylim',[0 600])
%%
print(['iE=',num2str(iE),' Event Count vs ShapeScoreNormalized_th zoom.tif'],'-dtiff');


%% (4.3 c) how [ShapeScoreNormalized_th] affect [ACC,precision,sensitivity], precision = PPV = TP/TP+FP, and sensitivity TPR=TP/TP+FN   
figure; hold on;
plot(scoreSorted, ACC,'linewidth',3,'color',[0 0 0]);
plot(scoreSorted, PPV,'linewidth',1.5,'color',[1 0 0]);
plot(scoreSorted, TPR,'linewidth',1.5,'color',[0 0 1]);

xlabel('\eta_{th}'); % xlabel('StrainScore Threshold');
ylabel('Rate');
set(gca,'fontsize',18,'xlim',[0 1],'ylim',[0,1.1]);
legend({'Accuracy (ACC)','Precision (PPV)','Sensitivity (TPR)'},'location','best');
%%
print(['iE=',num2str(iE),' Event Proportion vs ShapeScoreNormalized_th.tif'],'-dtiff');




%%
close all;


