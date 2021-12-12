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
% chenzhe, 2018-09-04
% Comment: looks like this code focus on the combining effect of the 2 classifiers. 
% Note that part (1) directly shows the result.  
% But actually, part (2) should be run first to find what values should be used in part (1) 

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

%% (1.0) Select an iE for analysis
iE = 4;

colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','trueTwinMap','scoreCF','cVolGrowthRatioMap','tProbMaxMap','mDistMap','sfMap');


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

mPhi = info(:,2);
mEta = info(:,3);
mDist = info(:,4);
mSF = info(:,5);
mPsize = info(:,6);
mPcnn = info(:,7);
T = table(response,mPhi,mEta,mDist,mSF,mPsize,mPcnn,  mDist.*mSF,  mDist.*mSF.*mPsize.*mPcnn,...
    'VariableNames',{'response','mPhi','mEta','mDist','SF','Psize','Pcnn','var1','var2'});
% save(['TwinS',num2str(iE),'.mat'],'response','predictors','T');

%% (1.1) Try fit a tree model
tc = fitctree(predictors,response,'MaxNumSplits',2);
view(tc,'mode','graph');


%% (1.2) Prepare Summarize how metrics change in the 2D-space of [x1-phi, x2-eta]
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

%% (1.3a) AND
TP = cumP;
FN = P-TP;
FP = cumN;
TN = N-FP;
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
ACC = (TP+TN)./(TP+FP+FN+TN);

ACC_max = max(ACC(:))
[ir,ic] = find(ACC == ACC_max, 1, 'first');
th = unique([x1(ir),x2(ic)],'rows')
TP_max = TP(sub2ind([nRC,nRC],ir,ic))
FP_max = FP(sub2ind([nRC,nRC],ir,ic))
FN_max = FN(sub2ind([nRC,nRC],ir,ic))
TN_max = TN(sub2ind([nRC,nRC],ir,ic))

[f,a,c] = myplot(X1,X2,ACC); title('ACC with AND relationship','fontweight','normal');
xlabel('\phi_{th}');ylabel('\eta_{th}');
set(a,'fontsize', 18);
caxis([0.87, 0.965]);
hold on; plot3(th(1),th(2),1,'xk','LineWidth',3,'MarkerSize',12);
%% Could compare naive algorithm

%%  (1.3b) use both score [phi_classifier (c=7 H=1)] AND [eta_classifier]
colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

strainScoreCF = 0.9503;    % 0.22, 0.17
shapeScoreCF = 0.0404;
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
            if (strainScore > strainScoreCF)&&(shapeScore > shapeScoreCF)
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else        % ground false
            if (strainScore > strainScoreCF)&&(shapeScore > shapeScoreCF)
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
confMap((strainScoreMapNew>strainScoreCF)&(shapeScoreMapNew>shapeScoreCF)) = confMap((strainScoreMapNew>strainScoreCF)&(shapeScoreMapNew>shapeScoreCF))+2;

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
ttl = ['(StrainScoreNormalized-th=',num2str(strainScoreCF),') AND (ShapeScoreNormalized-th=',num2str(shapeScoreCF),')'];
ttl_text = ['(\phi_{th}= ',num2str(strainScoreCF),') AND (\eta_{th}= ',num2str(shapeScoreCF),')'];
title(a,ttl,'fontweight','normal');
annotation(f,'textbox', [0.635 0.96 0.3 0.042], 'String',ttl_text, 'LineStyle','none', 'FontSize',24);
%% maximize plot and run this:
script_make_double_axis;
print([ttl,'.tif'],'-dtiff');


%% (1.4a)OR
TP = cumP(:,end)+cumP(end,:)-cumP;
FN = P-TP;
FP = cumN(:,end)+cumN(end,:)-cumN;
TN = N-FP;
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
ACC = (TP+TN)./(TP+FP+FN+TN);

ACC_max = max(ACC(:))
[ir,ic] = find(ACC == ACC_max, 1, 'first');
th = unique([x1(ir),x2(ic)],'rows')
TP_max = TP(sub2ind([nRC,nRC],ir,ic))
FP_max = FP(sub2ind([nRC,nRC],ir,ic))
FN_max = FN(sub2ind([nRC,nRC],ir,ic))
TN_max = TN(sub2ind([nRC,nRC],ir,ic))

[f,a,c] = myplot(X1,X2,ACC); title('ACC with OR relationship','fontweight','normal');
xlabel('\phi_{th}');ylabel('\eta_{th}');
set(a,'fontsize', 18);
caxis([0.87, 0.965]);
hold on; plot3(th(1),th(2),1,'xk','LineWidth',3,'MarkerSize',12);

%%  (1.4b) use both score [phi_classifier (c=7 H=1)] OR [eta_classifier]
colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

strainScoreCF = 0.9711;    % 0.22, 0.17
shapeScoreCF = 0.5370;
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
%%
close all;


