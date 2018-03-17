% confusion analysis for simScore criterion vs shapeScore criterion.


[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

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

%% (1)
colors = lines(7);
colorMap = [0 0 0; 0 0 1; colors(3,:); 1 0 0];
iE = 4;

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','simScoreMap','shapeScoreMap','trueTwinMap');


if strcmpi(sampleName,'WE43_T6_C1')
    switch iE
        case 5
            strainScoreCF = 0.230-0.5;
        case 4
            strainScoreCF = 0.215-0.5;
        case 3
            strainScoreCF = 0.217-0.5;
        case 2
            strainScoreCF = 0.291-0.5;
    end
end
% strainScoreCF = -0.5 + 0.285;   % s4, use default, or 0.284
shapeScoreCF = 0.1;
% previous criterion, 7Dis+(0.5-SF)< scoreCF
% i.e., new criterion, 7Dis-SF<socre-0.5

% select use strainScore or shapeScore
useStrainScore = 0;
useShapeScore = 0;
useBothScore = 1;


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





% plot map
map = zeros(size(exx));

% so that 3=TP(hit), 2=FP(false alarm), 1=FN(miss), 0=gb
map(trueTwinMap>0)=1;   % true twin

if useStrainScore
    % (1) simScore as criterion
    map(simScoreMap<strainScoreCF) = map(simScoreMap<strainScoreCF)+2;
elseif useShapeScore
    % (2) shapeScore as criterion
    map(shapeScoreMap>shapeScoreCF) = map(shapeScoreMap>shapeScoreCF)+2;
elseif useBothScore
    % (3) combine both as criterion
    map((simScoreMap<strainScoreCF)|(shapeScoreMap>shapeScoreCF)) = map((simScoreMap<strainScoreCF)|(shapeScoreMap>shapeScoreCF))+2;
end
% make only boundary as 0, others nan
map(map==0) = nan;
map(boundaryTFB==1) = 0;

[f,a,c] = myplot(map);
colormap(colorMap);
set(c,'limits',[0.75,3]);
c.Ticks = [1.125, 1.875, 2.625];
c.TickLabels={['FN(miss):',num2str(FN)],['FP(flase alarm):',num2str(FP)],['TP(hit):',num2str(TP)]}
set(a,'fontsize',18,'xticklabel',{''},'yticklabel',{''});
title(a,'');

%%

print('rename.tif','-dtiff');   % to parent folder

%%
for iE = 2:5

fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned','simScoreMap','shapeScoreMap','trueTwinMap');
twinSizePct(iE) = sum(trueTwinMap(:)>0)/sum(trueTwinMap(:)>=0)

end

