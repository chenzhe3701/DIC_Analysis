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
%
% chenzhe, 2018-05-09. Change confusion maps to use threshold value at
% optimized values.
%
% chenzhe, 2018-09-04.
% Based on script_2nd_paper_6_() and _7_(), find the threshold values at
% which the phi-classifier, and eta-classifier has the best performance, or
% their AND/OR relationship has the best performance.
%
% For more plots, at individually selected iEs, go to the script_ function

clear;clc;
addChenFunction;
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');    
end
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

f2 = '_';

printTF = 0;










%% This part is from code ()_script_2nd_paper_6_(), to find what is the [phi_th] ([eta_th]) for the best performance (high ACC value) of the [classifier] at [all] strain levels
%% (3) Confusion analysis using ACC metric, plot ACC curve vs criterion threshold
useACC = 1;
useF1 = 1-useACC;

%% (3.1) [ACC] vs. [phi_th] alone (C=7) , compare at [ALL] strain levels
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
if printTF
    print(['Event Proportion vs StrainScoreNormalized_th c_7.tif'],'-dtiff');
end

%% (3.2) [ACC] vs. [eta_th] alone, compare at [ALL] strain levels
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
if printTF
    print(['Event Proportion vs ShapeScoreNormalized_th.tif'],'-dtiff');
end





%% From the above, we can select some threshold values to observe/visualize
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
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','trueTwinMap','scoreCF','cVolGrowthRatioMap','tProbMaxMap','mDistMap','sfMap');
%% (1.1) Use [phi-classifier (c=7, H=1)], set (phi_th=0.9648)
strainScoreCF = 0.9648;    % 0.22, 0.17
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
%         shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
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
if printTF
    print([ttl,'.tif'],'-dtiff');
end


%% (1.3) use [eta-classifier], set (eta_th=0.5075)
strainScoreCF = 0.965;    % 0.22, 0.17
shapeScoreCF = 0.5075;

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
%         strainScore = stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H);
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
if printTF
    print([ttl,'.tif'],'-dtiff');
end


%% (1.4) use both score [phi-classifier (c=7 H=1)] OR [eta-classifier], set (phi_th=0.9648)OR(eta_th=0.5075)
strainScoreCF = 0.9648;    % 0.22, 0.17
shapeScoreCF = 0.5075;
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
if printTF
    print([ttl,'.tif'],'-dtiff');
end










%% This part is from ()_script_2nd_paper_7, to summarize what are the combinations of [phi_th] and [eta_th] when the two [classifiers] are used with [AND] or [OR] relationship.
%% (2) summarize at all strain levels, how the best ACC can be achieved with [And] and [Or] relationship, i.e., at what [phi_th, eta_th] combination
for iE = 2:5    
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap','scoreCF','cVolGrowthRatioMap','tProbMaxMap','mDistMap','sfMap');
    
    
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
            % [trueTwin, phi_normalized(C,H), eta_normalized, mDist, SF, P_size, P_cnn]
            info = [info; logical(stru(iS).trueTwin(iCluster)),...
                stru(iS).sf(iCluster)/(1+C*H) - C/(1+C*H)*min(H,stru(iS).dis(iCluster)) + (1+2*C*H)/2/(1+C*H),...
                (stru(iS).cVolGrowthRatio(iCluster,iE)+1)/2*stru(iS).tProbMax(iCluster),...
                stru(iS).dis(iCluster),...
                stru(iS).sf(iCluster),...
                stru(iS).cVolGrowthRatio(iCluster,iE),...
                stru(iS).tProbMax(iCluster),...
                stru(iS).dis(iCluster) * stru(iS).sf(iCluster),...
                stru(iS).dis(iCluster) * stru(iS).sf(iCluster) * stru(iS).cVolGrowthRatio(iCluster,iE) * stru(iS).tProbMax(iCluster)
                ];
        end
    end
    
    response = info(:,1);
    predictors = info(:,2:end);
    
    
    %% Prepare Summarize how metrics change in the 2D-space of [x1-phi, x2-eta]
    % chenzhe, comment on 2018-05-13. Looks like the 'sort' can be replaced by 'unique' to reduce redundancy. of course, other changes are required correspondingly, e.g., flpud, rep by numele, M++, etc  
    nRC = length(response);
    % ind1 is the row_number in info, which is unique between 1 and nR.
    % if we sort ind1 again, then rNum is the row_numbers in x1.
    [x1,ind1] = sort(info(:,2),'descend');  % phi-classifier intervals
    [~,rNum] = sort(ind1);
    [x2,ind2] = sort(info(:,3),'descend');  % eta-classifier
    [~,cNum] = sort(ind2);
    X1 = repmat(x1,1,nRC);
    X2 = repmat(x2',nRC,1);
    M = zeros(nRC);
    twins = [];
    for ii = 1:nRC
        if(1 == info(ii,1))
            M(rNum(ii),cNum(ii)) = 1;
        else
            M(rNum(ii),cNum(ii)) = -1;
        end
    end
    cumP = cumsum(cumsum(1==M,1),2);
    P = cumP(end);
    cumN = cumsum(cumsum(-1==M,1),2);
    N = cumN(end);
    
    %%  AND
    TP = cumP;
    FN = P-TP;
    FP = cumN;
    TN = N-FP;
    TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
    PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
    ACC = (TP+TN)./(TP+FP+FN+TN);
    
    ACC_max_and(iE) = max(ACC(:));
    [ir,ic] = find(ACC == ACC_max_and(iE), 1, 'first');
    th_and(iE,:) = unique([x1(ir),x2(ic)],'rows');
    TP_max_and(iE) = TP(sub2ind([nRC,nRC],ir,ic));
    FP_max_and(iE) = FP(sub2ind([nRC,nRC],ir,ic));
    FN_max_and(iE) = FN(sub2ind([nRC,nRC],ir,ic));
    TN_max_and(iE) = TN(sub2ind([nRC,nRC],ir,ic));
    TPR_max_and(iE) = TPR(sub2ind([nRC,nRC],ir,ic));
    PPV_max_and(iE) = PPV(sub2ind([nRC,nRC],ir,ic));
    
    %% OR
    TP = cumP(:,end)+cumP(end,:)-cumP;
    FN = P-TP;
    FP = cumN(:,end)+cumN(end,:)-cumN;
    TN = N-FP;
    TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
    PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
    ACC = (TP+TN)./(TP+FP+FN+TN);
    
    ACC_max_or(iE) = max(ACC(:));
    [ir,ic] = find(ACC == ACC_max_or(iE), 1, 'first');
    th_or(iE,:) = unique([x1(ir),x2(ic)],'rows');
    TP_max_or(iE) = TP(sub2ind([nRC,nRC],ir,ic));
    FP_max_or(iE) = FP(sub2ind([nRC,nRC],ir,ic));
    FN_max_or(iE) = FN(sub2ind([nRC,nRC],ir,ic));
    TN_max_or(iE) = TN(sub2ind([nRC,nRC],ir,ic));
    TPR_max_or(iE) = TPR(sub2ind([nRC,nRC],ir,ic));
    PPV_max_or(iE) = PPV(sub2ind([nRC,nRC],ir,ic));
    
end

t_and = table(th_and(iE_start:iE_stop,1),...
    th_and(iE_start:iE_stop,2),...
    ACC_max_and(iE_start:iE_stop)',...
    TP_max_and(iE_start:iE_stop)',...
    FP_max_and(iE_start:iE_stop)',...
    FN_max_and(iE_start:iE_stop)',...
    TN_max_and(iE_start:iE_stop)',...
    PPV_max_and(iE_start:iE_stop)',...
    TPR_max_and(iE_start:iE_stop)',...
    'VariableNames',{'Phi_th','Eta_th','Accuracy','TP','FP','FN','TN','Precision','Sensitivity'})

t_or = table(th_or(iE_start:iE_stop,1),...
    th_or(iE_start:iE_stop,2),...
    ACC_max_or(iE_start:iE_stop)',...
    TP_max_or(iE_start:iE_stop)',...
    FP_max_or(iE_start:iE_stop)',...
    FN_max_or(iE_start:iE_stop)',...
    TN_max_or(iE_start:iE_stop)',...
    PPV_max_or(iE_start:iE_stop)',...
    TPR_max_or(iE_start:iE_stop)',...
    'VariableNames',{'Phi_th','Eta_th','Accuracy','TP','FP','FN','TN','Precision','Sensitivity'})



