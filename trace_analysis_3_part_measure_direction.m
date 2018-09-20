
% chenzhe, 2018-09-09, around this time.
% 
% For each grain, we have some cluster number maps.  Maybe manually
% confirmed. But the true twin system number was not well identified,
% because we know that a pair of twin systems can exhibit similar surface
% strain.
%
% This code try to analyze the twin trace direction, to determine the
% active twin system. It should also be useful for slip trace analysis.

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
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
end
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

ID_target = 1144;

%% [data 1] Load cluster number data at all stops
struCell = cell(1,length(STOP)-1);

% [data 1] cluster data
for iE = iE_start:iE_stop
    
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    
    cNumMaps{iE} = clusterNumMap;
    cNumMaps_cleaned{iE} = clusterNumMapCleaned;
    
    struCell{iE} = stru;
end

%% [data 2] strain data. Convert into v7.3 for partial loading
clear strainFile;
for iE = iE_start:iE_stop
    strainFileName = [dicPath,'\',f2,STOP{iE+B}];
    disp(strainFileName);
    if ~exist([strainFileName,'_v73.mat'],'file')
        load(strainFileName);
        clear('exy_corrected');
        load(strainFileName,'exy_corrected');   % if 'exy_corrected' does not exist, this does not give error, rather, just warning.
        
        if exist('exy_corrected','var')&&(1==exy_corrected)
            disp('================= exy already corrected ! ========================');
            exy_corrected = 1;
        else
            disp('================= exy being corrected here ! =======================');
            exy = -exy;
            exy_corrected = 1;
        end
        % remove bad data points
        exx(sigma==-1) = nan;
        exy(sigma==-1) = nan;
        eyy(sigma==-1) = nan;
        qt_exx = quantile(exx(:),[0.0013,0.9987]); qt_exx(1)=min(-1,qt_exx(1)); qt_exx(2)=max(1,qt_exx(2));
        qt_exy = quantile(exy(:),[0.0013,0.9987]); qt_exy(1)=min(-1,qt_exy(1)); qt_exy(2)=max(1,qt_exy(2));
        qt_eyy = quantile(eyy(:),[0.0013,0.9987]); qt_eyy(1)=min(-1,qt_eyy(1)); qt_eyy(2)=max(1,qt_eyy(2));
        ind_outlier = (exx<qt_exx(1))|(exx>qt_exx(2))|(exy<qt_exy(1))|(exy>qt_exy(2))|(eyy<qt_eyy(1))|(eyy>qt_eyy(2));
        exx(ind_outlier) = nan;
        exy(ind_outlier) = nan;
        eyy(ind_outlier) = nan;
        
        outlier_removed = 1;
        save([strainFileName,'_v73.mat'],'outlier_removed','exy_corrected','-v7.3');
        
        myFile = matfile(strainFileName);
        myFields = who(myFile);
        for ii=1:length(myFields)
            save([strainFileName,'_v73.mat'],myFields{ii},'-append','-v7.3');
        end
    else
        disp('v7.3 file already exist');
    end
    strainFile{iE} = matfile([strainFileName,'_v73.mat']);
end

%% For each grain: () predict trace direction, () extract cluster map, () cleanup, () measure trace direction
clc;
iE = 4;
for iS = 1%:length(stru)
    %
    iS = find(arrayfun(@(x) x.gID == 378,stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace]
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
    traceSF = abs_schmid_factor(nss+1:nss+ntwin,2);
    
    % (2) Here, we want to do analysis on clustered maps, to decide the active slip/twin system.
    % For examples in twinning analysis, we previously performed cluster analysis and 'identified'/'confirmed' twin clusters.
    % More generally, we might need to first do a rough clustering based on strain map, then perform trace analysis, to decide which are the exist slip/twin systems.
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    
    % If need strain data (e,u,v,x,y), it should be re-imported for every new strain level
    needStrain = 0;
    if needStrain
        exx_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max); % strain of this region: grain + neighbor.
        exy_local = strainFile{iE}.exy(indR_min:indR_max, indC_min:indC_max);
        eyy_local = strainFile{iE}.eyy(indR_min:indR_max, indC_min:indC_max);
        sigma_local = strainFile{iE}.sigma(indR_min:indR_max, indC_min:indC_max);
    end
    
    clusterNumMapL = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    myplot(clusterNumMapL);
    
    % (3) Fill holes. Right now this is a little bit slow. May need to increase speed.
    clusterNumMapL = one_pass_fill(clusterNumMapL, 0.001);
    % myplot(clusterNumMapL);
    
    % (4) Decide if any is twin cluster?  Choose a/clusters.
    % --> Strain distribution might help to narrow down clusters that need analysis.
    % We mainly want to analyze those with high strains.
    clusterNumMapC = clusterNumMapL;    % for this cluster
    iC = 1;
    clusterNumMapC(clusterNumMapC~=iC) = 0;
    % myplot(clusterNumMapC);
    
    % (5) Then Do thinning/skeleton. The bwskel() function can perform some prunning at the same time.
    % Default no pruning. Because it uses 8-connectivity in bwskel_(), prunning sometimes makes analysis worse.
    thinTF = 1;    % This is mainly for debug, but looks like we should always consider thinning/skeleton.
    if thinTF
        clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
        % clusterNumMapT = double(thin(clusterNumMapC,inf));
        % clusterNumMapT = double(bwmorph((clusterNumMapC),'thin',inf) );
    end
    myplotm(clusterNumMapL, 'TF',clusterNumMapT, 'r', 1);
    caxis([-0.1, max(clusterNumMapL(:))+0.1]);
    
    % (6) Then do hough transform. H = intensity on [rhos, thetas] map
    [H,Theta,Rho] = hough(clusterNumMapT,'RhoResolution',1);    % H: d-1 is rho, d-2 is theta.
    
    % (7) Find peaks. Set a [neighborhood size], d_width = 5% map size, d_angle = 5 deg.
    maxNumPeaks = 32;
    peaks = houghpeaks(H, maxNumPeaks, 'Threshold', 0.3 * max(H(:)), 'NHoodSize',[round_odd(0.05*min(size(clusterNumMapL))),5] );
    peakAngles = Theta(peaks(:,2));
    peakStrength = H(sub2ind(size(H),peaks(:,1),peaks(:,2)));
    disp( table(peakAngles(:),peakStrength(:),'VariableNames',{'PeakAngles','PeakStrength'}) );
    
    % This is just to [illustrate] where the peak is in the hough space
    myplot(Theta,Rho,H);
    xlabel('theta'); ylabel('rho');
    axis normal;
    hold on;
    for k = 1:size(peaks,1)
        xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
        plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
    end
    
    % (8) Find lines. (A) Should keep gap small to prevent joining too many irrelavent parts. (2) MinLength should be decent
    % This is mainly to [illustrate] approximately, what are the peaks that were found ?
    lines = houghlines(clusterNumMapC, Theta, Rho, peaks, 'FillGap',5, 'MinLength', 0.1*min(size(clusterNumMapL)));
    myplotm(clusterNumMapL, 'r', 1);
    caxis([-0.1, max(clusterNumMapL(:))+0.1]);
    hold on;
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','g');
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','r');
        plot(xy(2,1),xy(2,2),'^','LineWidth',2,'Color','r');
    end
    
    % (9) Determine the active variant/ss by matching the peakAngles with traceND.
    % Use a [5 deg] threshold. Then, for those within valid angle range, make a score = SF/deltaAngle.
    activeVote = zeros(size(traceND));
    angleThreshold = 5;
    for ip = 1:length(peakAngles)
        dAngle = abs(traceND - peakAngles(ip));
        dAngle(dAngle > angleThreshold) = inf;
        score = traceSF./dAngle;
        % normalize
        if max(score)>0
            score = score/max(score);
        end
        activeVote = activeVote + score;
        
        % Another method to vote, but seems not as good.
        %        [val,ind] = min( abs(traceND - peakAngles(ip)) );
        %        if val < angleThreshold
        %           activeVote(ind) = activeVote(ind) + 1;
        %        end
    end
    activeSS = activeVote > 0.3*max(activeVote);    % Any one larger than 30% max vote is also selected
    disp(table(traceSF,traceND,activeVote,double(activeSS),'variableNames',{'Theoretical_SF','Theoretical_trace_ND', 'Vote', 'Active_SS_Identified'}));    
    
    
    
    % If there are more than one active slip/twin systems, we should do this:
    switch sum(activeSS)
        case 0
        case 1
        otherwise
    end
    
    % (9) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
    % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
    
    % (9.1) break skeleton into small branches
    skl = clusterNumMapT;
    branchPoints = bwmorph(skl, 'branchpoints');
    
    % can do this again to ensure to break skeleton, but not sure if its in general good ---------------
    heavyClean = 1;
    if heavyClean
        branchPoints = imdilate(branchPoints, ones(3));
    end
    branch = skl - branchPoints;
    
    % (9.2) assign an ID to each skeleton branch
    branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
    branchNumbered(~branch) = 0;
    % and get the unique ID of the branches
    uniqueBranchNum = unique(branchNumbered(:));
    uniqueBranchNum(uniqueBranchNum==0)=[];
    % [illustrate] skeleton branches
    myplotm(mod(branchNumbered,5)+logical(branchNumbered));
    
    % (9.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.  
    % Assign the ts/ss ID to the branches, which can be considered as grouped.  
    branchGrouped = zeros(size(branchNumbered));
    for ii = 1:length(uniqueBranchNum)
        model = fitlm(x_local(branchNumbered==uniqueBranchNum(ii)), y_local(branchNumbered==uniqueBranchNum(ii)));
        branchND = atand(-1/model.Coefficients.Estimate(2));
        dAngle = abs(traceND - branchND);
        dAngle(~activeSS) = inf;
        [~,ind] = min(dAngle);
        branchGrouped(branchNumbered == uniqueBranchNum(ii)) = nss + ind;
    end
    
    % (10) Grow each grouped branch into a a fragment with ID equals to active ss/ts.  
    [~,fragments] = city_block(branchGrouped);
    fragments(clusterNumMapC==0) = 0;
    % [illustrate] the fragments
    myplot(fragments, branch); caxis([18,24]);
    % ------------------ here, the one_pass_clean also need improvement. Basically, I want it not to allow to merge with NaN, but can into 0.  
%     if heavyClean
%         fragments = one_pass_fill(fragments, 0.01);
%         myplot(fragments, branch); caxis([18,24]);
%     end
    
    % (11 Next code maybe ?) 
    % Looks like at this stage, it can be considered as a twin system map that is more accurate than before.  Also, we have identified the active ts/ss. 
    % For each ts/ss frament/variant, analyze its position/connection with grain boundary, etc.
    
    
    
end

