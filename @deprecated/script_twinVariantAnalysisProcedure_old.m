% chenzhe, 2019-01-20
% use one grain as an example to show the process for data analysis including:
% (1) clean cluster number map  (cluster tracking? find gb cluster?)
% (2) cluster to twin by trace analysis
% (3) get twin gb intersection

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
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\*.mat','select the confirmed results where twin identification was based on trace dir and strain');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','gPhi1','gPhi','gPhi2');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','cityDistMap','ID','gID','gExx','gPhi1','gPhi','gPhi2');
end

% Make unique grain boundary map, and a list of the unique grain boundaries
[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:));
uniqueBoundaryList(uniqueBoundaryList==0)=[];

% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));


% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
useParallel = 1;

threshold = 1000;
useThreshold = 0;

% [data] strain data. Convert into v7.3 for partial loading
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

%% (0) load data: clusterNumberMapCell, stru-->struCell.  And show how to clean clusterNumMap
close all;
clusterNumberMapCell = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);

load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','trueTwinMapCell');     % not loaded 'twinMapCell','sfMapCell','tNote'

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMap;   % clusterNumMapCleaned;         % for this code use non-cleaned, and show how to clean
    % initialize/zero related fields
    for iS =1:length(struCell{iE})
        struCell{iE}(iS).cActiveSS = zeros(length(struCell{iE}(iS).cLabel), length(struCell{iE}(iS).tLabel));
    end
end


% This is clean cluster number map -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
iE_select = 4;
iC_select = 1;

for iE = iE_start:iE_stop
    
    iS = find(arrayfun(@(x) x.gID == 694,struCell{iE}));    % 378, 694  ----------------------------------------------------------------------------------------------------------------------------------------- 
    
    ID_current = gIDwithTrace(iS);
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    clusterNumMapLocal = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapLocal(ID_local~=ID_current) = 0;  % cluster number just this grain
    
    x_local = X(indR_min:indR_max, indC_min:indC_max);
    y_local = Y(indR_min:indR_max, indC_min:indC_max);
    boundaryTFLocal = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    exxLocal = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
    % (Figure-1)
    if ismember(iE,iE_select)
        exxLocal(ID_local~=ID_current)=nan;
        myplot(x_local, y_local, exxLocal, (boundaryTFLocal));
        title('\epsilon_x_x','fontweight','normal');
        set(gca,'xTick',[],'yTick',[],'fontsize',18);
        
        clusterNumMapLocal_a = clusterNumMapLocal;
        clusterNumMapLocal_a(ID_local~=ID_current)=nan;
        clusterNumMapLocal_a(clusterNumMapLocal_a==0)=nan;
        [f,a,c] = myplot(x_local, y_local, clusterNumMapLocal_a, (boundaryTFLocal));
        caxis([0 2]);
        set(c,'Ticks',[0 1 2]);
        title('Cluster ID','fontweight','normal');
        set(gca,'xTick',[],'yTick',[],'fontsize',18);
        colormap(parula(3));
        set(c,'limits',[2/3,2],'Ticks',[1,1+2/3],'TickLabels',{'1','2'});
    end
    % choices can be 0.0025 or 0.00025, whichever ok to keep more small and long traces. In real code = 0.001.  ------------------------------------------------- 
    clusterNumMapLocal = one_pass_fill_and_clean(clusterNumMapLocal, 0.00025);   
    if ismember(iE,iE_select)
        clusterNumMapLocal_a = clusterNumMapLocal;
        clusterNumMapLocal_a(ID_local~=ID_current)=nan;
        clusterNumMapLocal_a(clusterNumMapLocal_a==0)=nan;
        [f,a,c] = myplot(x_local, y_local, clusterNumMapLocal_a, (boundaryTFLocal));
        caxis([0 2]);
        set(c,'Ticks',[0 1 2]);
        title('Cluster ID','fontweight','normal');
        set(gca,'xTick',[],'yTick',[],'fontsize',18);
        colormap(parula(3));
        set(c,'limits',[2/3,2],'Ticks',[1,1+2/3],'TickLabels',{'1','2'});
    end
    
    % Because we have just cleaned clusterNumMap again, we need to update the clusterNumMapCell, which is different from the real code, which just load existing cleaned data
    % This is likely slow, but maybe enough for current purpose ?
    ind1 = find(ID==ID_current);
    ind2 = find(ID_local==ID_current);
    clusterNumberMapCell{iE}(ind1) = clusterNumMapLocal(ind2);
end


%% (1) cluster to twin by trace analysis
debugTF = 1;
iS = iS;
% clean old data to relabel.
for iE = iE_start:iE_stop
    struCell{iE}(iS).cActiveSS = zeros(size(struCell{iE}(iS).cActiveSS));
end

close all;
ID_current = struCell{iE}(iS).gID;
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

% initialize for each grain (iS)
for iE = iE_start:iE_stop
    twinMapLocal{iE} = zeros(size(ID_local));
    sfMapLocal{iE} = zeros(size(ID_local));
end
twinMapCell_cluster = []; % for cluster
sfMapCell_cluster = [];

% for each iE_entry (the entry point for analysis)
for iE_entry = iE_start:iE_stop
    % for each iC_entry
    for iC_entry = 1:length(struCell{iE_entry}(iS).cLabel)
        % We need to analyze this cluster [iC_outer] at the strain level [iE_outer], but this will need the information from the tracked [iE_list] and [iC_list].
        % So, first find the [iE_list, iC_list]
        [iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_entry, iC_entry);
        
        % Analyze all the linked iEs. Only need to analyze once: So, if iE_list(1)==iE_outer, it means it has not been analyzed before, then do [iE_list(ii),iC_list(ii)] pairs
        if iE_list(1) == iE_entry
            for iEC = 1:length(iE_list)
                %                 close all;
                iE = iE_list(iEC);
                iC = iC_list(iEC);
                
                %if(true) % This can control the actual selected strain level for anlaysis -------------------------------------------------------
                if (iE==iE_select)&&(iC==iC_select)&&(debugTF)  
                    ssAllowed = ones(ntwin,1);
                    script_label_twin_trace;    % -------------------------------------------------------------> this is the script to plot everything   
                    error('abc');
                end
                
                % run this for real analysis.  Below: each cell contains cells of tMap at an iEs
                ssAllowed = ones(ntwin,1);
                [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
                            struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,'twin',0 * debugTF, 0.3,0.3,ssAllowed, true);
                        
            end % end of iEC
        end
    end % end of iC_outer
end % end of iE_outer

%% for each strain level, update twinMapLocal{iE} with tMapCell
for iE = iE_start:iE_stop
    for jj = 1:size(twinMapCell_cluster,2)
        if (iE<=size(twinMapCell_cluster,1))&&(jj<=size(twinMapCell_cluster,2)) && (~isempty(twinMapCell_cluster{iE,jj}))
            twinMapLocal{iE} = twinMapLocal{iE} + twinMapCell_cluster{iE,jj};
            sfMapLocal{iE} = sfMapLocal{iE} + sfMapCell_cluster{iE,jj};
        end
    end
    if ~isempty(twinMapLocal{iE}) && (sum(twinMapLocal{iE}(:))>0)
        % myplot(twinMapLocal{iE}); title(['twinMap at strain: ',num2str(iE)]);   caxis([18 24]);
        % myplot(sfMapLocal{iE}); title(['sfMap at strain: ',num2str(iE)])
    end
    
    % Because we are just running example for one selected grain, we do not need to update the whole map of the sample, which is what we did in the real code.
    % But we might need to use the twinMap as the 'trueTwinMap' for analysis in the following section, finding twin/grain boundary intersection  ????
    ind1 = find(ID==ID_current);
    ind2 = find(ID_local==ID_current);
    trueTwinMapCell{iE}(ind1) = twinMapLocal{iE}(ind2);
end
disp(iS);


%% Update: I think now there is no need to demonstrate here, as it is already demonstrated in script_label_twin_trace().  The following used a different method, but no use, no change.
% %% (2) Then need to demonstrate how to find twin boundary intersection.  
% iE = iE_select;   %iS = find(arrayfun(@(x) x.gID == 694,struCell{iE}));    % 378, 694  --------------------------------
% close all;
% % Can initialize field
% struCell{iE}(iS).tGb = cell(1,ntwin);
% struCell{iE}(iS).tGbPts = cell(1,ntwin);
% for ii = 1:ntwin
%     struCell{iE}(iS).tGbPts{ii} = {[]};
% end
% 
% % If there are twins active: if sum(struCell{iE}(iS).cTrueTwin(:))>0
% if sum(struCell{iE}(iS).cTrueTwin(:))>0
%     
%     uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
%     uniqueBoundary_local = imdilate(uniqueBoundary_local,ones(3));
%     uniqueBoundary_local(ID_local~=ID_current)=0;
%     
%     trueTwinMapL = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
%     trueTwinMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
% 
%     trueTwinMapT = double( bwskel(imbinarize(trueTwinMapL),'MinBranchLength', 0 * round(min(size(trueTwinMapL))*0.2)) );
%     % new code to remove small branches from skeleton
%     [trueTwinMapT, branchPoints] = clean_skl(trueTwinMapT, 0.05*min(size(trueTwinMapL)) );
%     
%     skl_tn = trueTwinMapT .* trueTwinMapL;
%     % branchPoints = bwmorph(skl_tn, 'branchpoints');
%     branchPoints = imdilate(branchPoints, ones(3));
%     branch_tn = skl_tn;
%     branch_tn(branchPoints>0) = 0;
%     
%     branch_N = one_pass_label_8(branch_tn);    % here should use 8-connectivity to label
%     branch_N(~branch_tn) = 0;
%     branch_N(branchPoints) = 0;   % Note, this is new comparing to my old codes
% %     % Potential improvement, thin map for each individual twin, then combine the skeleton, but turned out to be a very bad method
% %     branch_N = zeros(size(trueTwinMapL));
% %     for it = nss+1 : nss+ntwin
% %         branch_N = branch_N + it * double( bwskel((trueTwinMapL==it),'MinBranchLength',0 * round(min(size(trueTwinMapL))*0.2)) );
% %     end    
% %     branch_N = one_pass_label_8(branch_N);
%     if debugTF
%         myplot(branch_N, boundaryTFLocal); caxis([0, max(branch_N(:))]);
%     end
%         
%     [~,fragments_N] = city_block(branch_N);
%     fragments_N(trueTwinMapL==0) = 0;
%     % Now at this point, [branch_tn] = short skeleton, numbered by the TS#.
%     % [branch_N] = short skeleton, numbered by sorting.
%     % [fragments_N] = twin fragments, numbered by sorting
%     if debugTF
%         myplot(fragments_N, boundaryTFLocal); caxis([0, max(branch_N(:))]);
%     end
%     
%     uniqueBranch_N = unique(branch_N(:));
%     uniqueBranch_N(uniqueBranch_N==0) = [];
%     % for each twin skeleton, need to have a minimum length along the trace direction
%     activeTS = sum(struCell{iE}(iS).cTrueTwin,1);
%     
%     dist_cr = 10;    % ------------------------------- This is a quite arbitrarily selected criterion  ---- -
%     npts_cr = 5;    % # of points within a certain distance to a gb ponint, then this frag can be considered as touching gb
%     
%     if debugTF
%         myplot(fragments_N,boundaryTFLocal); hold on;
%     end
%     for ii = 1:length(uniqueBranch_N)
%         bn = uniqueBranch_N(ii);
%         fragMap = fragments_N;          % choice-1: branch_N
%         fragMap(fragments_N~=bn) = 0;   % choice-1: branch_N
%         % check length use x_local and y_local
%         xrange = range(x_local(fragMap>0));
%         yrange = range(y_local(fragMap>0));
%         span = sqrt(xrange^2+yrange^2);
%         spanEnough = span > min(size(y_local,1),size(y_local,2))*(y_local(2)-y_local(1))*0.10;    % min(height,width) * scale * pct
%         
%         if spanEnough
%             tsNum = mode(trueTwinMapL(fragMap>0));    % most frequent element, is the ts number of this branch  % choice-1: branch_N
%             
%             iTwin = tsNum - nss;
%             
%             theta_target = traceND(iTwin);
%             [H,Theta,Rho] = hough(fragMap,'RhoResolution',1);
%             % block other peaks, if any
%             H(:,abs(Theta-theta_target)>2) = 0;
%             peaks = houghpeaks(H, 1);
%             peakRho = Rho(peaks(1));
%             peakTheta = Theta(peaks(2));
%             
%             lines = houghlines(uniqueBoundary_local, Theta, Rho, peaks, 'FillGap',999999);   % Note that the lines can be empty
%             
%             if ~isempty(lines)
%                 xy = [lines(1).point1; lines(1).point2];
% 
%                 % method-2: for each end point of houghline, calculate a distMap, and check if there is enough twin points within certain distance
%                 for jj = 1:2
%                     distMap = zeros(size(fragMap));
%                     distMap(xy(jj,2),xy(jj,1)) = 1;
%                     distMap = bwdist(distMap);
%                     % criterion for twin touch grain boundar: this fragment has > npts_cr ponints within distance dist_cr.
%                     if sum(sum( (distMap<dist_cr)&(fragMap>0) )) > npts_cr
%                         if debugTF
%                             plot3(xy(:,1),xy(:,2),[100;100],'LineWidth',1,'Color','g');
%                             % Plot beginnings and ends of lines
%                             plot3(xy(jj,1),xy(jj,2),[100;100],'o','LineWidth',2,'Color','r');
%                         end
%                         % This is a contact point. Record the position, and the twin system
%                         gbNum = uniqueBoundary_local(xy(jj,2),xy(jj,1));  % get gbNum
%                         ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
%                         
%                         % convert xy to absolute coordinate using x_local, y_local
%                         xcoord = x_local(xy(jj,2),xy(jj,1));
%                         ycoord = y_local(xy(jj,2),xy(jj,1));
%                         
%                         if isempty(ind)
%                             % append
%                             struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
%                             iGb = length(struCell{iE}(iS).tGb{iTwin});
%                             struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
%                         else
%                             % add
%                             iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
%                             struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
%     
%     %% do this for a specific fragment
%     if 0
%         myplot(fragments_N==38,boundaryTFLocal); hold on;
%         ii = 37;
%         
%         bn = uniqueBranch_N(ii);
%         fragMap = fragments_N;          % choice-1: branch_N
%         fragMap(fragments_N~=bn) = 0;   % choice-1: branch_N
%         % check length use x_local and y_local
%         xrange = range(x_local(fragMap>0));
%         yrange = range(y_local(fragMap>0));
%         span = sqrt(xrange^2+yrange^2);
%         spanEnough = span > min(size(y_local,1),size(y_local,2))*(y_local(2)-y_local(1))*0.10;    % min(height,width) * scale * pct
%         
%         if spanEnough
%             tsNum = mode(trueTwinMapL(fragMap>0));    % most frequent element, is the ts number of this branch  % choice-1: branch_N
%             
%             iTwin = tsNum - nss;
%             
%             theta_target = traceND(iTwin);
%             [H,Theta,Rho] = hough(fragMap,'RhoResolution',1);
%             % block other peaks, if any
%             H(:,abs(Theta-theta_target)>2) = 0;
%             peaks = houghpeaks(H, 1);
%             peakRho = Rho(peaks(1));
%             peakTheta = Theta(peaks(2));
%             
%             lines = houghlines(uniqueBoundary_local, Theta, Rho, peaks, 'FillGap',999999);   % Note that the lines can be empty
%             
%             if ~isempty(lines)
%                 xy = [lines(1).point1; lines(1).point2];
% 
%                 % method-2: for each end point of houghline, calculate a distMap, and check if there is enough twin points within certain distance
%                 for jj = 1:2
%                     distMap = zeros(size(fragMap));
%                     distMap(xy(jj,2),xy(jj,1)) = 1;
%                     distMap = bwdist(distMap);
%                     % criterion for twin touch grain boundar: this fragment has > npts_cr ponints within distance dist_cr.
%                     if sum(sum( (distMap<dist_cr)&(fragMap>0) )) > npts_cr
%                         if debugTF
%                             plot3(xy(:,1),xy(:,2),[100;100],'LineWidth',1,'Color','r');
%                             % Plot beginnings and ends of lines
%                             plot3(xy(jj,1),xy(jj,2),[100;100],'o','LineWidth',2,'Color','g');
%                         end
%                         % This is a contact point. Record the position, and the twin system
%                         gbNum = uniqueBoundary_local(xy(jj,2),xy(jj,1));  % get gbNum
%                         ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
%                         
%                         % convert xy to absolute coordinate using x_local, y_local
%                         xcoord = x_local(xy(jj,2),xy(jj,1));
%                         ycoord = y_local(xy(jj,2),xy(jj,1));
%                         
%                         if isempty(ind)
%                             % append
%                             struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
%                             iGb = length(struCell{iE}(iS).tGb{iTwin});
%                             struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
%                         else
%                             % add
%                             iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
%                             struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
%     
% end














