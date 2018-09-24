% clusterToTwin_modify

% chenzhe, 2018-09-19
%
% After running 3D_(), need some modification.
% For correction, select a grain, and correct all strain levels.

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

debugTF = 1;

%% [data] strain data. Convert into v7.3 for partial loading
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

%% (0) load data
cluster_number_maps = cell(1,length(STOP)-1);    % store all the clusterNumMap s, omit stop-0
cluster_number_maps_cleaned = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    cluster_number_maps{iE} = clusterNumMap;
    cluster_number_maps_cleaned{iE} = clusterNumMapCleaned;
    
end
% load struCell{}, twinMap{}
load('twinMaps.mat');

%% plot something of interest
close all;
iE_select = 2;
myplot(X,Y,twinMap{iE_select},boundaryTFB);caxis([18,24]);
myplot(X,Y,cluster_number_maps_cleaned{iE_select},boundaryTFB);

%% get ID from map, then
ids = find_ID_on_map(X,Y,ID,gcf,gca);

%% plot unit cell to check trace, one at a time
ind = find(gID==ids(1));
hcp_cell('euler',[gPhi1(ind),gPhi(ind),gPhi2(ind)], 'ss', 25:30, 'stress', [-1 0 0; 0 0 0; 0 0 0]);

%% (1) analyze
% get the local indices

stru = struCell{iE_start};
% for iS = 1 %1:length(stru)
% 181, 1350, 1390, 697, 193
iS = find(arrayfun(@(x) x.gID == 193,stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace], 
struCell{2}(iS).cVolGrowthRatio
%     iS = find(gIDwithTrace == 296); % for debugging.
% close all;
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

for iE = iE_start:iE_stop
    twinMapLocal{iE} = zeros(size(ID_local));
end

iLoop_iE = iE_start;
iLoop_iC = 1;
iLoop_iEC = 1;
cVolPctOld = 0;
cVolPctNotDecrease = 1;

%% This is the loop to run
close all;
% for iE_outer = iE_start:iE_stop
iE_outer = iLoop_iE;

% for each iE_outer, iC_outer, find the tracked iE_list, iC_list.

% for iC_outer = 1:length(struCell{iE_outer}(iS).cLabel)
iC_outer = iLoop_iC;

[iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_outer, iC_outer);

% Analyze all the linked iEs.  So, if iE_list(1)==iE_outer, it means it has not been analyzed before, then do [iE_list(ii),iC_list(ii)] pairs

iEC = iLoop_iEC;

% for iEC = 1:length(iE_list)

iE = iE_list(iEC);
iC = iC_list(iEC);
disp(['------------------------ [iE_outer, iC_outer, iE, iC] = [',num2str(iE_outer),', ',num2str(iC_outer),', ',num2str(iE),', ',num2str(iC),'] ------------------------']);

if iE_list(1) == iE_outer
    close all;
    clusterNumMapL = cluster_number_maps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    if debugTF
        myplot(clusterNumMapL);
    end
    
    % (3) Already filled holes.
    
    % (4) Decide if any is twin cluster?  Check if the strain is reasonable compared to theoretical twin strain.   
    % 
    % (4.1) Strain should be close, e.g., dist < 1.5 min_dist
    % (4.2) Effective strain should also be similar in magnitude. Otherwise, it is just closest, but not similar
    % (4.3) Maybe, select the clusters with higher strains, e.g. 1-out-of-2, 2-out-of-3, 2-out-of-4, 3-out-of-5, ...   
    % (4.4) Not active before
    cNum = struCell{iE}(iS).cLabel(iC);
    indClusterLocal = (clusterNumMapL==cNum);
    
    pdistCS = pdist2(struCell{iE}(iS).cCen(iC,:), struCell{iE}(iS).tStrain);
    minPdist = min(pdistCS);
    ok_1 = (pdistCS(:)<1.5*minPdist);
    
    eCluster = effective_strain_nBy3(struCell{iE}(iS).cCen(iC,:));
    eTwin = effective_strain_nBy3(struCell{iE}(iS).tStrain);
    ok_2 = (eCluster > 0.4*eTwin(:)) & (eCluster < 1.5*eTwin(:));
    
    cAllCluster = effective_strain_nBy3(struCell{iE}(iS).cCen);
    [sorted,rank_in_raw] = sort(cAllCluster);
    [tf,rank_in_sorted] = ismember(cAllCluster,sorted);
    rank_in_sorted_0_base = rank_in_sorted - 1; % convert to 0-based rank for easier comparison
    ok_3 = rank_in_sorted_0_base(iC) >= 0.3 * (length(struCell{iE}(iS).cLabel)-1);
    ok_3 = ones(size(ok_1)) * ok_3;
    
    %strainOKSS =  ok_1 & ok_2 & ok_3;
    
    disp(['eCluster = ', num2str(eCluster)]);
    disp(['all clusters strain = ',num2str(cAllCluster')]);
    disp(['higher strain cluster ?, ok_3 = ', num2str(ok_3(1))]);
    disp(table(pdistCS(:), double(ok_1(:)),eTwin(:), double(ok_2(:)),...
        'variableNames',{'pdistCS','ok_1','eTwin','ok_2'}));
    
    clusterNumMapC = clusterNumMapL;    % for this cluster.  -- Note that sometimes, the cluster was already cleaned to 0 size.
    clusterNumMapC(clusterNumMapC~=iC) = 0;
    cVolPct = sum(clusterNumMapC>0)/sum(clusterNumMapL>0);
    if cVolPct < cVolPctOld
        cVolPctNotDecrease = 0;
    else
        cVolPctNotDecrease = 1;
    end
    cVolPctOld = cVolPct;
    % myplot(clusterNumMapC);
    
    % (5) Then Do thinning/skeleton. The bwskel() function can perform some prunning at the same time.
    % Default no pruning. Because it uses 8-connectivity in bwskel_(), prunning sometimes makes analysis worse.
    thinTF = 1;    % This is mainly for debug, but looks like we should always consider thinning/skeleton.
    if thinTF
        clusterNumMapT = double( bwskel(imbinarize(clusterNumMapC),'MinBranchLength',0 * round(min(size(clusterNumMapC))*0.05)) );
        % clusterNumMapT = double(thin(clusterNumMapC,inf));
        % clusterNumMapT = double(bwmorph((clusterNumMapC),'thin',inf) );
    end
    if debugTF
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
    % disp( table(peakAngles(:),peakStrength(:),'VariableNames',{'PeakAngles','PeakStrength'}) );
    
    % This is just to [illustrate] where the peak is in the hough space
    if debugTF
        myplot(Theta,Rho,H);
        xlabel('theta'); ylabel('rho');
        axis normal;
        hold on;
        for k = 1:size(peaks,1)
            xy = [Theta(peaks(k,2)), Rho(peaks(k,1))];
            plot3(xy(1),xy(2),max(H(:)),'s','LineWidth',((maxNumPeaks+1-k)/maxNumPeaks)*4,'Color','k');
        end
    end
    
    % (8) Find lines. (A) Should keep gap small to prevent joining too many irrelavent parts. (2) MinLength should be decent
    % This is mainly to [illustrate] approximately, what are the peaks that were found ?
    if debugTF
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
    end
    
    % (9) Determine the active variant/ss by matching the peakAngles with traceND.
    % Use a [5 deg] threshold. Then, for those within valid angle range, make a score = SF/deltaAngle.
    traceVote = zeros(size(traceND));
    angleThreshold = 5;
    % [To prevent it is an all-zero map].
    if sum(peakStrength)>0
        for ip = 1:length(peakAngles)
            dAngle = abs(traceND - peakAngles(ip));
            dAngle(dAngle > angleThreshold) = inf;
            dAngle(dAngle < 1) = 1;
%             score = traceSF./dAngle;  % here we want to achieve that, for dAngle sasitfied, even if traceSF < 0, it still contributes
            score = logsig(transfer_to_logsig(traceSF, 0.2, 0.4, 0.9)) ./ dAngle;
            
            % normalize
            if max(score)>0
                score = score/max(score);
%             elseif min(score)<0
%                 % if there are traces match direction, but has negative SF
%                 score = (0.5-traceSF)./dAngle;
%                 score = score/max(score);
            end


            traceVote = traceVote + score;
        end
    end
    
    % [Need enough distinct peaks] The voted trace should be distinct. So if max(traceVote) < 0.3 * length(peakAngles), that means it's 'junk' vote
    enough_votes = max(traceVote) >= length(peakAngles) * 0.3;
    
    traceOKSS = (traceVote > 0.3*max(traceVote)) .* enough_votes;    % Any one larger than 30% max vote is also selected --> this need re-tunning
    
%     % [Additionally] If traces match super good, but clusterSize very small, then maybe it's ok. --------------------> This criterion need tunning.
%     % Where is the cluster 'small' criterion? .
    cVolPct = sum(clusterNumMapL(:)==iC)/sum(clusterNumMapL(:));
%     small_cluster_good_trace = zeros(size(traceOKSS));
%     [val, ind] = max(traceVote);
%     if val >= 0.75 * length(peakAngles)
%         small_cluster_good_trace(ind) = 1;
%     end
    
%     activeSS = strainOKSS & traceOKSS | small_cluster_good_trace;  % combine strainOK and traceOK
    
    % [10] Should also consider the activeSS from previous step, and combine
    if iEC > 1
        refActiveSS = struCell{iE_list(iEC-1)}(iS).cActiveSS(iC_list(iEC-1),:);
    else
        refActiveSS = zeros(ntwin,1);
    end
    ok_4 = ~refActiveSS(:);
    
    % ----------------------- criterion -----------------------------------------------------------------------------------------------------------------------------------------
    % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % traceMatch  
    % ok_1: strain close enough
    % ok_2: strain high enough
    % ok_3: it is the high strain cluster
    % refActiveSS: previously strain level active
    % ok_4: previous strain NOT active
    % cVolPctNotDecrease: cluster vol not decrease compared to previous strain level 
    activeSS = (traceOKSS & (ok_3 & ((ok_1&ok_2)|ok_4))) | refActiveSS(:);
    
    if debugTF
        disp(['# of peaks found: ', num2str(length(peakAngles))]);
        disp(['cluster vol pct: ', num2str(cVolPct)]);
        disp(['cluster vol not decrease: ', num2str(cVolPctNotDecrease)]);
        disp(table(traceSF,traceND,traceVote,double(traceOKSS),double(ok_1),double(ok_2),double(ok_3),double(ok_4),refActiveSS(:),double(activeSS),...
            'variableNames',{'traceSF','traceND','traceVote','traceOk','ok_1','ok_2','ok_3','ok_4','refActiveSS','activeSS'}));
    end
    
    % record the activeSS
    struCell{iE}(iS).cActiveSS(iC,:) = activeSS;
    
    % If there are more than one active slip/twin systems, we should seperate them:
    switch sum(activeSS)
        case 0
            fragments = zeros(size(clusterNumMapC));
        case 1
            ind = find(activeSS);
            fragments = (nss+ind) * ones(size(clusterNumMapC));
            fragments(clusterNumMapC==0) = 0;
            if debugTF
                myplot(fragments); caxis([18,24]);
            end
        otherwise
            % (11) Then break skeleton. Use broken skeleton as seed to grow, to fragment the cluster.
            % Here we have room to improve -- only the 'end' branches need to be seperated. Basically, we need traversal from end ponints.
            
            % (11.1) break skeleton into small branches
            skl = clusterNumMapT;
            branchPoints = bwmorph(skl, 'branchpoints');
            
            % can do this again to ensure to break skeleton, but not sure if its in general good ---------------
            heavyClean = 1;
            if heavyClean
                branchPoints = imdilate(branchPoints, ones(3));
            end
            branch = skl - branchPoints;
            
            % (11.2) assign an ID to each skeleton branch
            branchNumbered = one_pass_label_8(branch);    % here should use 8-connectivity to label
            branchNumbered(~branch) = 0;
            % and get the unique ID of the branches
            uniqueBranchNum = unique(branchNumbered(:));
            uniqueBranchNum(uniqueBranchNum==0)=[];
            % [illustrate] skeleton branches
            if debugTF
                myplotm(mod(branchNumbered,5)+logical(branchNumbered));
            end
            
            % (11.3) match each numbered skeleton branch to one of the active ts/ss, based on direction comparison.
            % Assign the ts/ss ID to the branches, which can be considered as grouped.
            branchGrouped = zeros(size(branchNumbered));
            % store the twin system r2 fit
            for itwin = 1:ntwin
                tR2{itwin} = 0;
            end
            for ib = 1:length(uniqueBranchNum)
                
                model = fitlm(x_local(branchNumbered==uniqueBranchNum(ib)), y_local(branchNumbered==uniqueBranchNum(ib)));
                
                branchND = atand(-1/model.Coefficients.Estimate(2));
                dAngle = abs(traceND - branchND);
                dAngle(~activeSS) = inf;
                [~,ind] = min(dAngle);
                branchGrouped(branchNumbered == uniqueBranchNum(ib)) = nss + ind;
                
                r2 = model.Rsquared.Ordinary;
                if isnan(r2)
                    r2 = 0;
                end
                tR2{ind} = [tR2{ind},r2];
            end
            for itwin = 1:ntwin
                struCell{iE}(iS).tR2(iC,itwin) = mean(tR2{itwin});
            end
            if debugTF
                struCell{iE}(iS).tR2(iC,:)
            end
            % tR2 only have number when multiple ts are fitted.
            
            
            % (12) Grow each grouped branch into a a fragment with ID equals to active ss/ts.
            [~,fragments] = city_block(branchGrouped);
            fragments(clusterNumMapC==0) = 0;
            % [illustrate] the fragments
            if debugTF
                myplot(fragments, branch); caxis([18,24]);
            end
    end
    %     if debugTF
    %         input(['finished iEC = ',num2str(iEC),', hit enter to continue']);
    %     end
    
    twinMapLocal{iE} = twinMapLocal{iE} + fragments;
    
    % end % end_For of iEC
    
end % end_If of if(iE_list(1)==iE_outer)

% update counter for next iteration
iLoop_iEC = iLoop_iEC + 1;
if iLoop_iEC > length(iE_list)
    iLoop_iEC = 1;
    iLoop_iC = iLoop_iC + 1;
    cVolPctOld = 0;
    cVolPctNotDecrease = 1;
end
if iLoop_iC > length(struCell{iE_outer}(iS).cLabel)
    iLoop_iEC = 1;
    iLoop_iC = 1;
    iLoop_iE = iLoop_iE + 1;
end
if iLoop_iE > iE_stop
    disp('finished iE,iC');
end

% end % end_For of iC_outer

% end % end_For iE_outer

% disp(iS);
% end % end of iS

%% update. First clean old map, then add new map.
for iE = iE_start:iE_stop
    toClean = twinMap{iE}(indR_min:indR_max, indC_min:indC_max);
    toClean(ID_local ~= ID_current) = 0;
    twinMap{iE}(indR_min:indR_max, indC_min:indC_max) = twinMap{iE}(indR_min:indR_max, indC_min:indC_max) - toClean + twinMapLocal{iE};
end


%%

timeStr = datestr(now,'yyyymmdd_HHMM');
save([timeStr,'_twinMaps.mat'],'twinMap','-v7.3');



