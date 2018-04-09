% confusion analysis for strainScore criterion vs shapeScore criterion.
%
% chenzhe, 2018-04-09 add note
% This code, based on manual labeled results, generate a 'confusion map'
% which plot clusters with different color to indicate the TP, FP, TN cases
% in classification. 


[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

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

%% (1) Generate a confusion map for selected iE. Label TP, FP, TN with different color. 
iE = 5;

colors = lines(7);
colorMap = [0 0 0; 0 0 1; colors(3,:); 1 0 0];

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap','scoreCF');

strainScoreCF = scoreCF;

% strainScoreCF = -0.5 + 0.285;   % s4, use default, or 0.284
shapeScoreCF = 0.1;
% previous criterion, 7Dis+(0.5-SF)< scoreCF
% i.e., new criterion, 7Dis-SF<socre-0.5

% select use strainScore or shapeScore
useStrainScore = 1;
useShapeScore = 0;
useBothScore = 0;


TP = 0;    % hit
FP = 0;    % false alarm
FN = 0;    % missing
TN = 0;    % correct rejection
if useStrainScore
    for iS=1:length(stru)
        for iCluster=1:length(stru(iS).cLabel)
            strainScore = 7*stru(iS).dis(iCluster)-stru(iS).sf(iCluster);
            shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
            if stru(iS).trueTwin(iCluster)>0    % if ground true
                if strainScore < strainScoreCF
                    TP = TP + 1;
                else
                    FN = FN + 1;
                    disp(iS)
                end
            else        % ground false
                if strainScore < strainScoreCF
                    FP = FP + 1;
                else
                    TN = TN + 1;
                end
            end
            
        end
    end
elseif useShapeScore
    for iS=1:length(stru)
        for iCluster=1:length(stru(iS).cLabel)
            strainScore = 7*stru(iS).dis(iCluster)-stru(iS).sf(iCluster);
            shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
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
elseif useBothScore
    for iS=1:length(stru)
        for iCluster=1:length(stru(iS).cLabel)
            strainScore = 7*stru(iS).dis(iCluster)-stru(iS).sf(iCluster);
            shapeScore = stru(iS).cvInc(iCluster)*stru(iS).tProbMax(iCluster);
            if stru(iS).trueTwin(iCluster)>0    % if ground true
                if (strainScore < strainScoreCF)||(shapeScore > shapeScoreCF)
                    TP = TP + 1;
                else
                    FN = FN + 1;
                end
            else        % ground false
                if (strainScore < strainScoreCF)||(shapeScore > shapeScoreCF)
                    FP = FP + 1;
                else
                    TN = TN + 1;
                end
            end
        end
    end
end



% plot 'confusion map'
confMap = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
confMap(trueTwinMap>0)=1;   % true twin

% add by 2. So if it's TP, it will be 1+2=3. If it's FP, it will be 0+2=2.
if useStrainScore
    % (1) strainScore as criterion
    confMap(strainScoreMap<strainScoreCF) = confMap(strainScoreMap<strainScoreCF)+2;
elseif useShapeScore
    % (2) shapeScore as criterion
    confMap(shapeScoreMap>shapeScoreCF) = confMap(shapeScoreMap>shapeScoreCF)+2;
elseif useBothScore
    % (3) combine both as criterion
    confMap((strainScoreMap<strainScoreCF)|(shapeScoreMap>shapeScoreCF)) = confMap((strainScoreMap<strainScoreCF)|(shapeScoreMap>shapeScoreCF))+2;
end
% make only boundary as 0, others nan
confMap(confMap==0) = nan;
confMap(boundaryTFB==1) = 0;

[f,a,c] = myplot(confMap);
colormap(colorMap);
set(c,'limits',[0.75,3]);
c.Ticks = [1.125, 1.875, 2.625];
c.TickLabels={['FN(miss):',num2str(FN)],['FP(flase alarm):',num2str(FP)],['TP(hit):',num2str(TP)]}
set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
title(a,'');

%%

print('rename.tif','-dtiff');   % to parent folder

%% (2) summarize the twin area/vol fraction at strain levels 2-5
for iE = 2:5    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','strainScoreMap','shapeScoreMap','trueTwinMap');
    twinSizePct(iE) = sum(trueTwinMap(:)>0)/sum(trueTwinMap(:)>=0);    
end

figure;
plot(strainPauses(2:5),twinSizePct(2:5)*100,'-or','linewidth',1.5);
set(gca,'xdir','reverse','fontsize',18);
xlabel('Global Uniaxial Strain, mm/mm');
ylabel('Twin Area Percent, %');

