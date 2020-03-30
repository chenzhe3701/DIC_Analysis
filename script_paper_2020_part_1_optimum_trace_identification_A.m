% This is for generating data

%% This is same as beginning of 'trace_analysis_3D_clusterToTwin_by_trace_strain.m'
clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2');
end
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

debugTF = 0;

%% (0) load data
% store all the clusterNumMap s, omit stop-0
% cluster_number_maps = cell(1,length(STOP)-1);    
clusterNumberMapCell = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
    twinMapCell{iE} = zeros(size(clusterNumMapCleaned));
    sfMapCell{iE} = zeros(size(clusterNumMapCleaned));
    cToGbDistMapCell{iE} = zeros(size(clusterNumMapCleaned));
    % initialize/zero related fields
    for iS =1:length(stru)
        stru(iS).cActiveSS = zeros(length(stru(iS).cLabel), length(stru(iS).tLabel));
    end
    % try to remove some fields, if needed
    try
        stru = rmfield(stru,{'tR2'});
    end
    struCell{iE} = stru;
end
%% we need to get ground truth data

[truthFile, truthPath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\*.mat','select the results for twin-grain boundary intersection');
load(fullfile(truthPath,truthFile),'struCell');

%%
iE_start = 2;
stru = struCell{iE_start};
[~, ~, nss, ntwin, ~] = define_SS(sampleMaterial,'twin');
t = [];
variableNames = {'ID','iE','iC','cTrueTwin','activePre','activePost',...
    'dStrain','rStrain','dEffStrain','rEffStrain','cStrainRank','cStrainRankMax',...
    'tSF','d1','d2','d3','d4','d5','nPeaks','traceDir','nd1','nd2','nd3','nd4','nd5'};
T = cell2table(cell(0,length(variableNames)));
T.Properties.VariableNames = variableNames;

for iS = 1:length(stru)
    %     iS = find(arrayfun(@(x) x.gID == 378,stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace]
    %     iS = find(gIDwithTrace == 296); % for debugging.
    close all;
    ID_current = stru(iS).gID;
    
    % (1) Calculate theoretical trace direction.
    ind_euler = find(gID==ID_current);
    euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
    if (1==eulerAligned)
        % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
    else
        % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
    end
    [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
    traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
    mult_factor = ones(size(traceDir));
    mult_factor(traceDir<0) = 1;
    mult_factor(traceDir>=0) = -1;
    traceND = traceDir + 90*mult_factor;    % convert traceDir to traceND
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    
    for iE_entry = iE_start:iE_stop
        % for each iC_outer
        for iC_entry = 1:length(struCell{iE_entry}(iS).cLabel)
            % We need to analyze this cluster [iC_outer] at the strain level [iE_outer], but this will need the information from the tracked [iE_list] and [iC_list].
            % So, first find the [iE_list, iC_list]
            [iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_entry, iC_entry);
            
            % Analyze all the linked iEs.  So, if iE_list(1)==iE_outer, it means it has not been analyzed before, then do [iE_list(ii),iC_list(ii)] pairs
            if iE_list(1) == iE_entry
                
                % This is the start of analyzing a new linked clusters
                for iEC = 1:length(iE_list)
                    close all;
                    iE = iE_list(iEC);
                    iC = iC_list(iEC);
                    
                    clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
                    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
                    
                    cTrueTwin = struCell{iE}(iS).cTrueTwin(iC,:);   % 1x6
                    cTrueTwin = cTrueTwin(:);
                    tSF = struCell{iE}(iS).tSF;   % 1x6
                    tSF = tSF(:);   % 6x1
                    tStrain = struCell{iE}(iS).tStrain;   % 6x3
                    tEffStrain = effective_strain_nBy3(tStrain); % 6x1
                    cStrain = struCell{iE}(iS).cCen(iC,:);    % 1x3
                    cEffStrain = effective_strain_nBy3(cStrain);   % 1x1 
                    
                    dStrain = pdist2(cStrain, tStrain);
                    dStrain = dStrain(:);
                    rStrain = dStrain/min(dStrain);
                    
                    dEffStrain = abs(cEffStrain - tEffStrain);
                    rEffStrain = cEffStrain./tEffStrain;
                    
                    cEffStrainAll = effective_strain_nBy3(struCell{iE}(iS).cCen);
                    [sorted,rank_in_raw] = sort(cEffStrainAll);
                    [tf,rank_in_sorted] = ismember(cEffStrainAll,sorted);
                    rank_in_sorted_0_base = rank_in_sorted - 1; % convert to 0-based rank for easier comparison
                    cStrainRank = rank_in_sorted_0_base(iC);
                    cStrainRankMax = max(rank_in_sorted_0_base);
                    
                    % Known to be active in later step. (If not every ss is allowed, then it must be a check coming back from a later step)
                    % When coming back, we might be able to relax requirement for ok_1 and ok_2
                    if iEC < length(iE_list)
                        active_post = reshape(struCell{iE_list(iEC+1)}(iS).cTrueTwin(iC_list(iEC+1),:), ntwin, 1);
                    else
                        active_post = zeros(ntwin,1);
                    end
                    %strainOKSS =  ok_3 & ((ok_1 & ok_2)|(active_post));
                    
                    % Know to be active in a previous step
                    if iEC > 1
                        active_pre = struCell{iE_list(iEC-1)}(iS).cActiveSS(iC_list(iEC-1),:);
                        active_pre = active_pre(:);
                    else
                        active_pre = zeros(ntwin,1);
                    end
                    
                    clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
                    clusterNumMapC(clusterNumMapC~=iC) = 0;
                    
                    
                    % (5) Then Do thinning/skeleton. The bwskel() function can perform some prunning at the same time.
                    clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
                    [clusterNumMapT, branchPoints] = clean_skl(clusterNumMapT, round(min(size(clusterNumMapC))*0.05));
                    if debugTF==1
                        myplotm(clusterNumMapL, 'TF',clusterNumMapT, 'r', 1);
                        caxis([-0.1, max(clusterNumMapL(:))+0.1]);
                    end
                    
                    % (6) Then do hough transform. H = intensity on [rhos, thetas] map
                    [H,Theta,Rho] = hough(clusterNumMapT,'RhoResolution',1);    % H: d-1 is rho, d-2 is theta.
                    
                    % (7) Find peaks. Set a [neighborhood size], d_width = 5% map size, d_angle = 5 deg.
                    maxNumPeaks = 32;
                    peaks = houghpeaks(H, maxNumPeaks, 'Threshold', 0.3 * max(H(:)), 'NHoodSize',[round_odd(0.05*min(size(clusterNumMapL))),5] );
                    peakAngles = Theta(peaks(:,2));
                    peakStrength = H(sub2ind(size(H),peaks(:,1),peaks(:,2)));
                    if debugTF >= 1
                        disp( table(peakAngles(:),peakStrength(:),'VariableNames',{'PeakAngles','PeakStrength'}) );
                    end
                    nPeaks = size(peaks,1);
                    
                    % This is just to [illustrate] where the peak is in the hough space
                    if debugTF==1
                        myplot(Theta,Rho,H);
                        xlabel('theta'); ylabel('rho');
                        axis normal;
                        hold on;
                        for k = 1:size(peaks,1)
                            xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
                            plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
                        end
                    end
                    
                    % Count #peaks within [0:5] degrees range? ==> Note: not cummulative.   
                    angDiffMat = zeros(6,5);   % angle diff matrix = nTwins x degree range  
                    if sum(peakStrength)>0
                        for ip = 1:length(peakAngles)
                            % Can preset a threshold, but maybe just record the SF rather than using a threshold
                            dAngle = min(abs(traceND - peakAngles(ip)), abs(traceND + 180*(-1)*sign(traceND) - peakAngles(ip)));    % ==> Important for calculating angle diff.
                            dAngleCeil = ceil(dAngle);
                            
                            % if multiple traces within angle of 5,  select the one with higher SF
                            if sum(dAngleCeil<=5)>1
                                ind = (dAngleCeil<=5);
                                maxSF = max(tSF(ind));
                                ind = (tSF==maxSF);
                                dAngleCeil(~ind) = dAngleCeil(~ind) + 5; %angleThreshold;
                            end
                            % end of treating multiple traces with close direction   
                            
                            for itwin = 1:6
                                if dAngleCeil(itwin)<=5
                                    angDiffMat(itwin,dAngleCeil(itwin)) = angDiffMat(itwin,dAngleCeil(itwin))+1;
                                end
                            end
                        end
                    else
                        nPeaks = nPeaks;
                        angDiffMat = angDiffMat;
                    end
                    
                    % count for each twin system, how many twin system have theoretical trace dir within [0: 1,2,3,4,5] degrees ==> cummulative 
                    twinAngleDiffMat = zeros(6,5);
                    for itwin = 1:6
                       tdir = traceDir(:);
                       tdir(itwin) = [];
                       angDiff = min(abs(traceDir(itwin)-tdir), abs(traceDir(itwin) + 180*(-1)*sign(traceDir(itwin)) - tdir));
                       dAngleCeil = ceil(angDiff);
                       for dAngle = 1:5
                          twinAngleDiffMat(itwin,dAngle) = sum(dAngleCeil<=dAngle); 
                       end
                    end
                    
                    % This summarize the strain criterion
                    t_local = [repmat([ID_current, iE, iC],6,1), cTrueTwin, active_pre, active_post, ...
                        dStrain, rStrain, dEffStrain, rEffStrain, repmat(cStrainRank,6,1), repmat(cStrainRankMax,6,1), ...
                        tSF(:), angDiffMat, repmat(nPeaks,6,1), traceDir, twinAngleDiffMat];
                    t = [t;t_local];
                    
                end
            end
        end
    end
    disp(['iS=',num2str(iS)]);
end
T=array2table(t);
T.Properties.VariableNames=variableNames;
%%
save('temp_results/optimum_identification.mat','T')