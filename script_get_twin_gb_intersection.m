% Note this is the code I try to determine the point at which twin touches
% the grain boundareis, etc.... 2018-11-07.
%
% chenzhe, 2018-10-20
%
% After we got an idea of suitable SF

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

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\','select the confirmed results where twin identification was based on trace dir and strain');

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

debugTF = 0;

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

%% (0) load data, using SF threshold values to assign active twin system, and make maps
% Load cluster number maps (cleaned).
clusterNumberMapCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');


%%
try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'cityDistMap','gNeighbors');
end

% Create a list that makes a note of the [iE, grainID, tSF, truceTwin
% Possibly need these data:
% (1) gNeighbors is up to date, because it can be updated in the code '_1C_EBSDtoSEM_after_GB_adjustment.m'.
% (2) cityDistMap (from gb), this was made to be saved in code '_3C_find_gb_cluster.m'

% Make unique grain boundary map, and a list of the unique grain boundaries
[~, boundaryID, neighborID, tripleTF, tripleID] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:));
uniqueBoundaryList(uniqueBoundaryList==0)=[];

%%
% The relatively easy and strainght-forward task, similar to paper
% For each twin, we want to find its intersecting gb, then look at the:
% (1) misorientation, (2) displacement gradient component in neighbor
plotTF = 0;
[~, ~, ~, ntwin, ~] = define_SS(sampleMaterial,'twin');

tb_gbNum = [];
tb_iE = [];
tb_gNum = [];
tb_tsNum = [];
tb_pts = [];


% ids = find_ID_on_map(X,Y,ID,gcf,gca);
% iS_target = find(arrayfun(@(x) x.gID == ids(1),struCell{2}))
% plotTF = 1;

for iE = iE_start:iE_stop
%%
for iS = 1:length(struCell{iE})     % e.g., iS=1253, ID_current=1390
    % Can initialize field
    struCell{iE}(iS).tGb = cell(1,ntwin);
    struCell{iE}(iS).tGbPts = cell(1,ntwin);
    for ii = 1:ntwin
        struCell{iE}(iS).tGbPts{ii} = {[]};
    end
    
    ID_current = struCell{iE}(iS).gID
    % If there are twins active: if sum(struCell{iE}(iS).cTrueTwin(:))>0
    
    % Calculate theoretical trace direction.
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
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    % Make it one data point wider on each side
    indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
    indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
    indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
    indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
    uniqueBoundary_local = imdilate(uniqueBoundary_local,ones(3));
    uniqueBoundary_local(ID_local~=ID_current)=0;
    
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    
    
    trueTwinMapL = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    trueTwinMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    
    trueTwinMapT = double( bwskel(imbinarize(trueTwinMapL),'MinBranchLength',0 * round(min(size(trueTwinMapL))*0.05)) );
    skl_tn = trueTwinMapT .* trueTwinMapL;
    branchPoints = bwmorph(skl_tn, 'branchpoints');
    branchPoints = imdilate(branchPoints, ones(3));
    branch_tn = skl_tn;
    branch_tn(branchPoints>0) = 0;
    
    branch_N = one_pass_label_8(branch_tn);    % here should use 8-connectivity to label
    branch_N(~branch_tn) = 0;
    branch_N(branchPoints) = 0;   % Note, this is new comparing to my old codes
    [~,fragments_N] = city_block(branch_N);
    fragments_N(trueTwinMapL==0) = 0;
    % Now at this point, [branch_tn] = short skeleton, numbered by the TS#.
    % [branch_N] = short skeleton, numbered by sorting.
    % [fragments_N] = twin fragments, numbered by sorting
    
    uniqueBranch_N = unique(branch_N(:));
    uniqueBranch_N(uniqueBranch_N==0) = [];
    % for each twin skeleton, need to have a minimum length along the trace direction
    activeTS = sum(struCell{iE}(iS).cTrueTwin,1);
    
    dist_cr = 10;    % ------------------------------- This is a quite arbitrarily selected criterion  ---- -
    npts_cr = 5;    % # of points within a certain distance to a gb ponint, then this frag can be considered as touching gb
    
    if plotTF
        myplotm(fragments_N, 'r', 1); hold on;
    end
    for ii = 1:length(uniqueBranch_N)
        bn = uniqueBranch_N(ii);
        fragMap = fragments_N;          % choice-1: branch_N
        fragMap(fragments_N~=bn) = 0;   % choice-1: branch_N
        % check length use x_local and y_local
        xrange = range(x_local(fragMap>0));
        yrange = range(y_local(fragMap>0));
        span = sqrt(xrange^2+yrange^2);
        spanEnough = span > min(size(y_local,1),size(y_local,2))*(y_local(2)-y_local(1))*0.10;    % min(height,width) * scale * pct
        
        if spanEnough
            tsNum = mode(trueTwinMapL(fragMap>0));    % most frequent element, is the ts number of this branch  % choice-1: branch_N
            
            iTwin = tsNum - nss;
            
            theta_target = traceND(iTwin);
            [H,Theta,Rho] = hough(fragMap,'RhoResolution',1);
            % block other peaks, if any
            H(:,abs(Theta-theta_target)>2) = 0;
            peaks = houghpeaks(H, 1);
            peakRho = Rho(peaks(1));
            peakTheta = Theta(peaks(2));
            
            lines = houghlines(uniqueBoundary_local, Theta, Rho, peaks, 'FillGap',999999);   % Note that the lines can be empty
            
            if ~isempty(lines)
                xy = [lines(1).point1; lines(1).point2];
                
                % method-1: calculate the distance map from the fragment. Then check the two end points of the houghline.
                %                distMap = bwdist(fragMap);
                %                for jj = 1:2
                %                   if distMap(xy(jj,2),xy(jj,1)) < dist_cr
                %                       % This is a contact point. Record the position, and the twin system
                %
                %                   end
                %                end
                
                % method-2: for each end point of houghline, calculate a distMap, and check if there is enough twin points within certain distance
                for jj = 1:2
                    distMap = zeros(size(fragMap));
                    distMap(xy(jj,2),xy(jj,1)) = 1;
                    distMap = bwdist(distMap);
                    if sum(sum( (distMap<dist_cr)&(fragMap>0) )) > npts_cr
                        if plotTF
                            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','g');
                            % Plot beginnings and ends of lines
                            plot(xy(jj,1),xy(jj,2),'x','LineWidth',2,'Color','r');
                        end
                        % This is a contact point. Record the position, and the twin system
                        gbNum = uniqueBoundary_local(xy(jj,2),xy(jj,1));  % get gbNum
                        ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                        if isempty(ind)
                            % append
                            struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                            iGb = length(struCell{iE}(iS).tGb{iTwin});
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  xy(jj,:);  % assign point coord (1x2 vector) to the cell value
                        else
                            % add
                            iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [struCell{iE}(iS).tGbPts{iTwin}{iGb};xy(jj,:)];    % append the point coord (1x2 vector) as new rows
                        end
                    end
                end
                
            end
        end
    end
    
    
    % So we need to correlate twin with its touching uniqueGB (therefore neighbor grain)
    % It's better to calculate a point of intersection, by modeling the gb and
    % twin both as a line, or line + eliptical axis
    % But later, we may also need to determine if two twins coincide at a gb.
    %
    % Algorithm:
    % For each uniqueGB, find the related grain(1) Take cluster number map, skeleton,
    % Growth to match size threshold, eliminate gb cluster
    % Fit to match direction criterion
    
end

%% 
tboundary = zeros(size(uniqueBoundary));

for iS = 1:length(struCell{iE})
    
    ID_current = struCell{iE}(iS).gID
    
    activeTS = find(sum(struCell{iE}(iS).cTrueTwin, 1)> 0);
    for ii = 1:length(activeTS)
        iTwin = activeTS(ii);
        tsNum = iTwin + nss;
        
        % find unique grain id, grain ID, etc
        GBs = struCell{iE}(iS).tGb{iTwin};
        Pts = struCell{iE}(iS).tGbPts{iTwin};
        for jj = 1:length(GBs)
            ind = find((tb_gbNum==GBs(jj))&(tb_iE==iE)&(tb_gNum==ID_current)&(tb_tsNum==tsNum));
            
            if isempty(ind)
                ind = size(tb_gbNum,1)+1
                tb_gbNum(ind,1) = GBs(jj);
                tb_iE(ind,1) = iE;
                tb_gNum(ind,1) = ID_current;
                tb_tsNum(ind,1) = tsNum;
                tb_pts{ind} = Pts{jj};
                tboundary(uniqueBoundary == GBs(jj)) = 1;
            else
                disp('error');
            end
            
        end
    end
end

%% plot
% myplot(strainFile{iE}.exx);limit = caxis; close;
% myplot(X,Y,strainFile{iE}.exx + grow_boundary(grow_boundary(logical(uniqueBoundary))),grow_boundary(grow_boundary(tboundary)));
% caxis(limit);
% myplot(X,Y,trueTwinMapCell{iE}+15*grow_boundary(grow_boundary(logical(uniqueBoundary))) ,grow_boundary(grow_boundary(tboundary)));


tBoundaryCell{iE} = tboundary;

end
%%
timeStr = datestr(now,'yyyymmdd_HHMM');
save([timeStr,'_twin_at_boundary_result.mat'], 'tb_gbNum', 'tb_iE', 'tb_gNum', 'tb_tsNum', 'tb_pts', 'tBoundaryCell','struCell','-v7.3');


%% Assume the previous analysis is correct.
% (1) calculate the average misorientation of twinned boundary vs non twinned boundary

figure;hold on;
for iE = 2:5
    gbTwinned_iE = tb_gbNum(tb_iE==iE);
    gbMiso = zeros(length(uniqueBoundaryList),1);
    ind_twinned = zeros(length(uniqueBoundaryList),1);
    
    for ii = 1:length(uniqueBoundaryList)
        gbNum = uniqueBoundaryList(ii);
        g1 = floor(gbNum/10000);
        g2 = mod(gbNum,10000);
        ind_euler = find(gID==g1);
        euler_1 = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
        ind_euler = find(gID==g2);
        euler_2 = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
        gbMiso(ii) = calculate_misorientation_hcp(euler_1,euler_2);
        if ismember(gbNum,gbTwinned_iE)
            ind_twinned(ii) = 1;
        end
    end
    
    gbMisoTwinned = gbMiso(ind_twinned==1);
    gbMisoNonTwinned = gbMiso;
    gbMisoNonTwinned(ind_twinned==1) = [];
    
    edges = 0:5:95;
    
    % figure;histogram(gbMisoTwinned,edges);
    % title(['With twin, iE=',num2str(iE)]);
    % xlabel('gb misorientation');
    % ylabel('counts');
    %
    % figure;histogram(gbMisoNonTwinned,edges);
    % title(['Without twin, iE=',num2str(iE)]);
    % xlabel('gb misorientation');
    % ylabel('counts');
    
    plot(edges(1:end-1), histcounts(gbMisoTwinned,edges)./(histcounts(gbMisoTwinned,edges)+histcounts(gbMisoNonTwinned,edges)), '-o')

end

title(['# grains twinned over # grains total']);
legend({'iE=2','iE=3','iE=4','iE=5'})
xlabel('gb misorientation');
ylabel(' ');

%% 








