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

%% (1.0) Select an iE for analysis
iE = 4;

colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

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
mPhi = info(:,2);
mEta = info(:,3);
mDist = info(:,4);
mSF = info(:,5);
mPsize = info(:,6);
mPcnn = info(:,7);
T = table(response,mPhi,mEta,mDist,mSF,mPsize,mPcnn,mDist.*mSF, mDist.*mSF.*mPsize.*mPcnn,...
    'VariableNames',{'response','mPhi','mEta','mDist','SF','Psize','Pcnn','var1','var2'});
% save(['TwinS',num2str(iE),'.mat'],'response','predictors','T');

%% (1.1) Fit some model
tc = fitctree(predictors,response,'MaxNumSplits',3);
view(tc,'mode','graph');






%% (1.2) Summarize how metrics change in the 2D-space of [x1-phi, x2-eta]
nRC = length(x1);
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

%% AND
TP = cumP;
FN = P-TP;
FP = cumN;
TN = N-FP;
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
ACC = (TP+TN)./(TP+FP+FN+TN);
myplot(X1,X2,ACC); title('AND');

ACC_max = max(ACC(:))
[ir,ic] = find(ACC == ACC_max);
th = unique([x1(ir),x2(ic)],'rows')
TP_max = TP(sub2ind([nRC,nRC],ir,ic))
FP_max = FP(sub2ind([nRC,nRC],ir,ic))
FN_max = FN(sub2ind([nRC,nRC],ir,ic))
TN_max = TN(sub2ind([nRC,nRC],ir,ic))
%% OR
TP = cumP(:,end)+cumP(end,:)-cumP;
FN = P-TP;
FP = cumN(:,end)+cumN(end,:)-cumN;
TN = N-FP;
TPR = TP./(TP+FN);  % true positive rate, sensitivity, recall. (How many of the positives are identified)
PPV = TP./(TP+FP);  % precision, positive predictive value. (How many identified are the real positives).
ACC = (TP+TN)./(TP+FP+FN+TN);
myplot(X1,X2,ACC);  title('OR')

ACC_max = max(ACC(:))
[ir,ic] = find(ACC == ACC_max);
th = unique([x1(ir),x2(ic)],'rows')
TP_max = TP(sub2ind([nRC,nRC],ir,ic))
FP_max = FP(sub2ind([nRC,nRC],ir,ic))
FN_max = FN(sub2ind([nRC,nRC],ir,ic))
TN_max = TN(sub2ind([nRC,nRC],ir,ic))



%% naive algorithm








%%  (1.4) use both score [StrainScoreNormalized (c=7 H=1)] OR [ShapeScore_New], set (phi_th=0.9496)AND(eta_th=0.0395), (0.9696)OR(0.5238)
colors = lines(7);
colorMap = [0 0 0; 1 1 1; 0 0 1; colors(5,:); 1 0 0];

strainScoreCF = 0.9696;    % 0.22, 0.17
shapeScoreCF = 0.5238;
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
%%






%% (0.5) plot strainScoreMapNormalized, adjust color so that high vs low scores are distinct
C = 7; H = 1;
strainScoreMapNew = sfMap/(1+C*H) -C/(1+C*H)*min(H,mDistMap) + (1+2*C*H)/2/(1+C*H);
[f,a,c] = myplot(X,Y,strainScoreMapNew,boundaryTFB,3);
% set(gca,'fontsize',24,'XTick',[],'YTick',[])
title('StrainScoreNormalized','fontweight','normal');
caxis([0.9 1]);
% colormap(parula(8));
%% maximize plot and run this:
script_make_double_axis;
title(c,'\phi','fontsize',24);
print('strain score Normalized iE=4.tif','-dtiff');



%% (0.6) plot ShapeScoreMapNormalized, adjust color so that high vs low scores are distinct
shapeScoreMapNew = (cVolGrowthRatioMap+1)/2.*tProbMaxMap;
[f,a,c] = myplot(X,Y,shapeScoreMapNew,boundaryTFB,3);
% set(gca,'fontsize',24,'XTick',[],'YTick',[])
title('ShapeScoreNormalized','fontweight','normal');
caxis([0 1]);
% colormap(parula(8));
%% maximize plot and run this:
script_make_double_axis;
title(c,'\eta','fontsize',24);
print('shape score Normalized iE=4.tif','-dtiff');



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




%% (3) Confusion analysis using ACC metric, plot ACC curve vs criterion threshold

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


