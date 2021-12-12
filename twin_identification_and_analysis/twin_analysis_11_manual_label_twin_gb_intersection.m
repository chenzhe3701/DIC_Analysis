% chenzhe, 2019-07-17
% It is likely useful to manually label twin-gb intersection.  We need some ground truth data anyway.


clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1\SEM Data\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('D:\p\m\DIC_Analysis\setting_for_real_samples\WE43_T6_C1_setting.mat','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1\Analysis_2021_09','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\WE43_T6_C1\Analysis_2021_09\WE43_T6_C1_3_twin_result_with_variant.mat','select the results where twin identification was based on trace dir and strain');

[twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\WE43_T6_C1\Analysis_2021_09\WE43_T6_C1_4_twin_result_with_gb_intersection.mat','select the results for twin-grain boundary intersection (auto analyzed)');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
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
%%
clear strainFile;
for iE = iE_start-1:iE_stop
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
for iE = []%iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
end

% Load from the pre-labeled results: twinMapCell, sfMapCell, struCell.  (cToGbDistMapCell is omitted, as will no longer be used in this code)
% load(fullfile(confirmedLabelPath,confirmedLabelFile),'struCell','twinMapCell','trueTwinMapCell','sfMapCell','tNote');
load(fullfile(confirmedLabelPath,confirmedLabelFile),'trueTwinMapCell');

% load previous twin_gb interaction result, for reference.
load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile));
% load('D:\WE43_T6_C1\Analysis_2021_09\possibly useful data\20191102_2147_twin_at_boundary_manual_result.mat','struCell');    % 2021-10-04 debug, for struCell_ref

%% [The following is for manual label]
struRef = struCell; % save the old struCell as ref.



%% Initialize the related fields
% Select an iE.
for iE = iE_start:iE_stop
    iS = 1;
    % iS = 162;
    % ID_current = 190;
    for iS = 1:length(struCell{iE})
        %% If juse run the following, it can clear the fields. -------------
        struCell{iE}(iS).tGb = cell(1,6);
        struCell{iE}(iS).tGbPts = cell(1,6);
        struCell{iE}(iS).tGbNormal = cell(1,6);
        struCell{iE}(iS).tGbPtsiE = cell(1,6);
        for ii = 1:6
            struCell{iE}(iS).tGbPts{ii} = {[]};
            struCell{iE}(iS).tGbNormal{ii} = {[]};
            struCell{iE}(iS).tGbPtsiE{ii} = {[]};
        end
        %% Can use this to initiate strain levels > 2, if higher level is based on lower level
        if (false)&&(iE>2)
            struCell{iE}(iS).tGb = struCell{iE-1}(iS).tGb;
            struCell{iE}(iS).tGbPts = struCell{iE-1}(iS).tGbPts;
            struCell{iE}(iS).tGbNormal = struCell{iE-1}(iS).tGbNormal;
            struCell{iE}(iS).tGbPtsiE = struCell{iE-1}(iS).tGbPtsiE;
        end
    end
end

%% Find twinned_grain_list
twinned_grain_list = [];
% Select an iE.
for iE = iE_start:iE_stop
    twinned_grain_list{iE} = [];
    for iS = 1:length(struCell{iE})
        if sum(struCell{iE}(iS).cTrueTwin(:))>0
            twinned_grain_list{iE} = [twinned_grain_list{iE}; struCell{iE}(iS).gID];
        end
    end
end

%% [for modification only] If modify earlier strain level (e.g., 2), also modify that grain at later strain level (copy 2->3, then 3->4, etc)
% % start to modify from this iE
% start_modify_from_iE = 5;
% for ie = start_modify_from_iE:iE_stop
% %     a0 = load(fullfile(saveDataPath, ['temp_results\for_recover_iS_',num2str(iS)]),'struCell');
% %     a0 = a0.struCell{ie-1}(iS);
%     a0 = struCell{ie-1}(iS);
%     % copy
%     struCell{ie}(iS).tGb = a0.tGb;
%     struCell{ie}(iS).tGbPts = a0.tGbPts;
%     struCell{ie}(iS).tGbNormal = a0.tGbNormal;
%     struCell{ie}(iS).tGbPtsiE = a0.tGbPtsiE;
% end


%% Select an iS to start
iS = 1;
%%
lines_total = 0;

iE = iE_start;
continueTF = true;
while (continueTF)&&(iS<length(struCell{iE}))
    close all;
    plottedRelatedInfo = false;
    
    for iE = iE_start:iE_stop                
        disp(['iE=',num2str(iE),', iS=',num2str(iS),', ID=',num2str(struCell{iE}(iS).gID)]);
        
        % [relabel 20191101 code] in this code, may need to initialize the 4 fields based on previous iE result. 
%         if iE>iE_start
%             struCell{iE}(iS).tGb = struCell{iE-1}(iS).tGb;
%             struCell{iE}(iS).tGbPts = struCell{iE-1}(iS).tGbPts;
%             struCell{iE}(iS).tGbNormal = struCell{iE-1}(iS).tGbNormal;
%             struCell{iE}(iS).tGbPtsiE = struCell{iE-1}(iS).tGbPtsiE;
%         end
            
        % find grains that have been determined to be twinned up to this iE
        gns = [];
        for ii = iE_start:iE
            gns = [gns;twinned_grain_list{ii}];
        end
        activeTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
        
        if ismember(struCell{iE}(iS).gID, gns)            
            %% close all;
            ID_current = struCell{iE}(iS).gID;
            ind = find(gID==ID_current);
            
            euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
            if (1==eulerAligned)
                % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
                [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
            else
                % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
                [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
            end
            [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
            traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
            
            nNeighbors = gNNeighbors(ind);
            ID_neighbors = gNeighbors(ind, 1:nNeighbors);
            
            ind_local = ismember(ID, [ID_current, ID_neighbors]); %ismember(ID, [ID_current,ID_neighbor]);
            
            % Make it one data point wider on each side
            indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
            indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
            indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
            indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
            X_local = X(indR_min:indR_max, indC_min:indC_max);
            Y_local = Y(indR_min:indR_max, indC_min:indC_max);
            uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
            boundaryTF_local = boundaryTF(indR_min:indR_max, indC_min:indC_max);
            trueTwinMapLocal = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
            e_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
            
            % plot only once, twinMap and eMap at related iEs
            if ~plottedRelatedInfo
                for ii = iE_stop:-1:iE+1
                    myplot(X_local, Y_local, trueTwinMapCell{ii}(indR_min:indR_max, indC_min:indC_max), grow_boundary(boundaryTF_local));
                    caxis([18 24]);
                    title(['iE=',num2str(ii)],'fontweight','normal');
                    disableDefaultInteractivity(gca);
                    
                    myplot(X_local, Y_local, strainFile{ii}.exx(indR_min:indR_max, indC_min:indC_max), grow_boundary(boundaryTF_local));
                    title(['iE=',num2str(ii)],'fontweight','normal');
                    for ID_target = ID_neighbors
                        if ismember(ID_target, unique(ID_local(:)))
                            local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, [1,2,3,19:24], gca);
                        end
                    end
                    disableDefaultInteractivity(gca);
                end
                plottedRelatedInfo = true;
            end
            
            [hF_t,a0,~] = myplot(X_local, Y_local, trueTwinMapLocal, boundaryTF_local);
            caxis([0 6]);
            label_map_with_ID(X_local, Y_local, ID_local, hF_t, ID_current);
            disableDefaultInteractivity(a0);
            
            [hF_e,aA,~] = myplot(X_local, Y_local, e_local, grow_boundary(boundaryTF_local));
            title(['iE=',num2str(iE)],'fontweight','normal');
            disableDefaultInteractivity(aA);
            hold on;
            active_ts_to_label = 18 + find(sum(struCell{iE}(iS).cTrueTwin,1));
            local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_current, active_ts_to_label, gca);
            %         for ID_target = ID_neighbors
            %             local_fun_label_map_with_ID_and_trace(X_local,Y_local,ID_local,ID_target, [1,2,3,19:24], gca);
            %         end
            
            
            %%
%             % --> find existing lines in this region
%             for iTwin = 1:6
%                 tGbNormals = struCell{iE}(iS).tGbNormal{iTwin};
%                 for iGb = 1:length(tGbNormals)
%                     if ~isempty(tGbNormals{iGb})
%                         for iLine = 1:length(tGbNormals{iGb})
%                             drawline(aA,'Position',tGbNormals{iGb}{iLine});
%                         end
%                     end
%                 end
%             end
            
            % [relabel 20191101 code] Input new trace, but here can get from reference/ground truth struCell 
            % An improvement could be, compare with lines at previous iE. only draw new lines. 
            % This is because manual label also have errors: sometimes a line with iE at current is actually an old line.  
            lines_previous_iE = {[0 0; 0 0]};
            ilines = 0;
            for ie = iE_start:iE
                for itwin = 1:6
                    for igb = 1:length(struCell{ie}(iS).tGb{itwin})
                        for ipt = 1:length(struCell{ie}(iS).tGbPtsiE{itwin}{igb})
                            if struCell{ie}(iS).tGbPtsiE{itwin}{igb}(ipt) <= iE
                                line_t = struCell{ie}(iS).tGbNormal{itwin}{igb}{ipt};    % tentative line
                                if all(cellfun(@(x) sum(abs(x(:)-line_t(:)))>50, lines_previous_iE))    % tentative line significantly different from any previous line
                                    ilines = ilines + 1;
                                    lines_previous_iE{ilines} = line_t;
                                    lines_total = lines_total + 1;
                                    drawline(aA,'Position', lines_previous_iE{ilines});     % draw it
                                end
                            end
                        end
                    end
                end
            end

%             lines_to_draw = [];
%             ilines = 0;
%             for itwin = 1:6
%                 for igb = 1:length(struCell{iE}(iS).tGb{itwin})
%                     for ipt = 1:length(struCell{iE}(iS).tGbPtsiE{itwin}{igb})
%                         if struCell{iE}(iS).tGbPtsiE{itwin}{igb}(ipt) == iE
%                             % if not found in previous iE
%                             line_t = struCell{iE}(iS).tGbNormal{itwin}{igb}{ipt};    % tentative line
%                             if all(cellfun(@(x) sum(abs(x(:)-line_t(:)))>50, lines_previous_iE))    % tentative line significantly different from any previous line
%                                 ilines = ilines + 1;
%                                 lines_to_draw{ilines} = line_t;
%                                 lines_total = lines_total + 1
%                             end
%                         end
%                     end
%                 end
%             end
%             
            needToAdd = questdlg('add trace?','select answer','Yes','No','Cancel','No');
            
%             % [relabel 20191101 code] assign needToAdd based on if there were lines to draw
%             if isempty(lines_to_draw)
%                 needToAdd = 'No';
%             else
%                 needToAdd = 'Yes';
%                 nlines = length(lines_to_draw);
%                 ilines = 1;
%             end
            
            switch needToAdd
                case 'Yes'
                    nextTrace = true;
                    while nextTrace
                        % can use a button to start drawing a line
                        handleDrawline = drawline(aA,'Color','r');
                        posAdd = customWait(handleDrawline);
                        
                        answer = questdlg('add,redraw,or next grain?','select operation',...
                            'NextTrace','Redraw','AcceptAndNext_iE',...
                            'AcceptAndNext_iE');

%                         % [relabel 20191101 code] assign the lines to draw to posAdd, and asign value to answer  
%                         posAdd = lines_to_draw{ilines};
%                         if ilines<nlines
%                             answer = 'NextTrace';
%                             ilines = ilines + 1;
%                         else
%                             answer = 'AcceptAndNext_iE';
%                         end
                        

                        % find coordinate/maybe indices
                        pt1 = posAdd(1,:);
                        pt2 = posAdd(2,:);
                        x = posAdd(:,1);
                        y = posAdd(:,2);
                        
                        % find intersection between the line drawn and uniqueGrainBoundary.
                        [pts,inds,indR,indC] = grids_covered_by_line(X_local,Y_local,pt1,pt2);
                        gbNum = uniqueBoundary_local(inds);
                        gbNum = gbNum(gbNum>0);
                        if isempty(gbNum)
                            % We need to make it a little bit thicker to ensure intersection
                            disp('Shift gb normal to make intersection');
                            posAdd(1,1) = posAdd(1,1) + X_local(1,2)-X_local(1,1);
                            posAdd(2,1) = posAdd(2,1) + X_local(1,2)-X_local(1,1);
                            pt1 = posAdd(1,:);
                            pt2 = posAdd(2,:);
                            x = posAdd(:,1);
                            y = posAdd(:,2);
                            [pts,inds,indR,indC] = grids_covered_by_line(X_local,Y_local,pt1,pt2);
                            gbNum = uniqueBoundary_local(inds);
                            gbNum = gbNum(gbNum>0);
                        end
                        % sometimes, the line covers two unique grain boundaries
                        igb = 1;
                        ngb = length(gbNum);
                        while igb<=ngb
                            gbn = gbNum(igb);
                            if (floor(gbn/10000)~=ID_current)&&(mod(gbn,10000)~=ID_current)
                                igb = igb + 1;
                            else
                                gbNum = gbn;
                                break;
                            end
%                             error("redraw, line covers more than one unique grain boundaries");
                        end
                        % [relabel 20191101 code]. Note, modify, but don't remember if this will work for the original code.  
                        % Today there is an example of  gbNum = [12691259; 12691259; 12951269]; @ iE=2, iS=1168, ID_curren = 1295.  
                        if igb>ngb
                            error("redraw, line covers more than one unique grain boundaries");
                        end
                                       
                        ind_of_ind = find(uniqueBoundary_local(inds)==gbNum, 1, 'first');
                        ind = inds(ind_of_ind);
                        xcoord = X_local(ind);
                        ycoord = Y_local(ind);
                        
                        switch answer
                            case 'NextTrace'
                                % look at direction of line drawn, determine which twin system it is
                                direction = atand((y(2)-y(1))/(x(2)-x(1)));
%                                 [~,iTwin] = min(abs(direction-traceDir));  

                                % [relabel 20191101 code] chenzhe, 2019-11-01
                                % a potential problem is that, if two systems have similar direction, 
                                % but the one with similar direction as the labeled was not determined as active, 
                                % it can lead to the labeling of a ts that we did not want to label.   
                                %
                                % Method: first, assign 'inf' to the dirDiff of the inactive twin systems.  Need to consider angle range. 
                                dirDiff = abs(direction-traceDir);
                                dirDiff(dirDiff>90) = 180 - dirDiff(dirDiff>90);
                                dirDiff(~activeTS) = inf;
                                [ang_min,iTwin] = min(dirDiff);
                                
                                % make sure the drawn line is similar to the direction of an active trace, otherwise do not add   
                                if ang_min < 20
                                    iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                                    if isempty(iGb)
                                        % append
                                        struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                                        iGb = length(struCell{iE}(iS).tGb{iTwin});
                                        struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                                        struCell{iE}(iS).tGbNormal{iTwin}{iGb} = {[pt1;pt2]};
                                        struCell{iE}(iS).tGbPtsiE{iTwin}{iGb} = iE;
                                    else
                                        % add
                                        iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                                        struCell{iE}(iS).tGbPts{iTwin}{iGb} = [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                                        struCell{iE}(iS).tGbNormal{iTwin}{iGb} = [struCell{iE}(iS).tGbNormal{iTwin}{iGb}; {[pt1;pt2]}];
                                        struCell{iE}(iS).tGbPtsiE{iTwin}{iGb} = [struCell{iE}(iS).tGbPtsiE{iTwin}{iGb}; iE];
                                    end
                                else
                                    disp('does not match active twin system, ignore');
                                    delete(handleDrawline);
                                end
                            case 'Redraw'
                                delete(handleDrawline);
                            case 'AcceptAndNext_iE'
                                % look at direction of line drawn, determine which twin system it is
                                direction = atand((y(2)-y(1))/(x(2)-x(1)));
%                                 [~,iTwin] = min(abs(direction-traceDir));
                                
                                % [relabel 20191101 code] chenzhe, 2019-11-01
                                % a potential problem is that, if two systems have similar direction, 
                                % but the one with similar direction as the labeled was not determined as active, 
                                % it can lead to the labeling of a ts that we did not want to label.   
                                %
                                % Method: first, assign 'inf' to the dirDiff of the inactive twin systems 
                                dirDiff = abs(direction-traceDir);
                                dirDiff(dirDiff>90) = 180 - dirDiff(dirDiff>90);
                                dirDiff(~activeTS) = inf;
                                [ang_min,iTwin] = min(dirDiff);
                                
                                if ang_min < 20
                                    iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                                    if isempty(iGb)
                                        % append
                                        struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                                        iGb = length(struCell{iE}(iS).tGb{iTwin});
                                        struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                                        struCell{iE}(iS).tGbNormal{iTwin}{iGb} = {[pt1;pt2]};
                                        struCell{iE}(iS).tGbPtsiE{iTwin}{iGb} = iE;
                                    else
                                        % add
                                        iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                                        struCell{iE}(iS).tGbPts{iTwin}{iGb} = [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                                        struCell{iE}(iS).tGbNormal{iTwin}{iGb} = [struCell{iE}(iS).tGbNormal{iTwin}{iGb}; {[pt1;pt2]}];
                                        struCell{iE}(iS).tGbPtsiE{iTwin}{iGb} = [struCell{iE}(iS).tGbPtsiE{iTwin}{iGb}; iE];
                                    end
                                else
                                    disp('does not match active twin system, ignore');
                                    delete(handleDrawline);
                                end
                                close(hF_e);
                                close(hF_t);
                                nextTrace = false;
                                
                                % if accept addition and then next iE, need to first copy fields at this iE to later iEs   
                                % if iE>=3
                                for ii = iE+1:iE_stop
                                    struCell{ii}(iS).tGb = struCell{iE}(iS).tGb;
                                    struCell{ii}(iS).tGbPts = struCell{iE}(iS).tGbPts;
                                    struCell{ii}(iS).tGbNormal = struCell{iE}(iS).tGbNormal;
                                    struCell{ii}(iS).tGbPtsiE = struCell{iE}(iS).tGbPtsiE;
                                end
                                % end   % [relabel 20191101 code] don't remember why used iE>3. Was that for debugging?  Disable it for now today.    
                        end
                    end
                    
                case 'No'
                    % if go to next iE, and this iE has labels not contained in next iE. need to first copy fields at this iE to later iEs  
                    % use 'isempty' because if just look and check, we should not copy 
                    if (iE<iE_stop)&&(isempty(cell2mat(struCell{iE+1}(iS).tGb)))
                        for ii = iE+1:iE_stop
                            struCell{ii}(iS).tGb = struCell{iE}(iS).tGb;
                            struCell{ii}(iS).tGbPts = struCell{iE}(iS).tGbPts;
                            struCell{ii}(iS).tGbNormal = struCell{iE}(iS).tGbNormal;
                            struCell{ii}(iS).tGbPtsiE = struCell{iE}(iS).tGbPtsiE;
                        end
                    end
                    
                    close(hF_e);
                    close(hF_t);
                case 'Cancel'
                    iS = iS - 1;    % when cancel, reduce by 1 first
                    %                 close(handleFigA);
                    %                 close(handleFig0);
                    continueTF = false;
                    break;
            end
            
        end
        
        
    end % end of iE=2:5
    
    iS = iS + 1;
end

%% it's better to save periodically
try
    save(fullfile(saveDataPath, ['for_recover_iS_',num2str(iS)]),'struCell','iE','iS','-append');
catch
    save(fullfile(saveDataPath, ['for_recover_iS_',num2str(iS)]),'struCell','iE','iS');
end

%% Record data at all iEs
iE_current = iE;    % record

% initiate
tb_gbNum = [];
tb_iE = [];
tb_gNum = [];
tb_tsNum = [];
tb_pts = [];
tb_iE_inGrain = []; % first iE level appear in this grain
tb_iE_atBoundary = [];  % first iE level appear at this boundary

tboundary = zeros(size(uniqueBoundary));

for iE = iE_start:iE_current
    % summarize
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID
        activeTS = find(sum(struCell{iE}(iS).cTrueTwin, 1)> 0);     % but maybe haven't intersected a boundary yet
        for ii = 1:length(activeTS)
            iTwin = activeTS(ii);
            tsNum = iTwin + 18;
            
            % find unique grain id, grain ID, etc
            GBs = struCell{iE}(iS).tGb{iTwin};
            Pts = struCell{iE}(iS).tGbPts{iTwin};
            iEs = struCell{iE}(iS).tGbPtsiE{iTwin};
            
            iE_inGrain = min(cell2mat(struCell{iE}(iS).tGbPtsiE{iTwin}(:)));
            
            for iGb = 1:length(GBs)
                iE_atBoundary = min(struCell{iE}(iS).tGbPtsiE{iTwin}{iGb});
                
                ind = find((tb_gbNum==GBs(iGb))&(tb_iE==iE)&(tb_gNum==ID_current)&(tb_tsNum==tsNum));
                if isempty(ind)
                    ind = size(tb_gbNum,1)+1
                    tb_gbNum(ind,1) = GBs(iGb);
                    tb_iE(ind,1) = iE;
                    tb_gNum(ind,1) = ID_current;
                    tb_tsNum(ind,1) = tsNum;
                    tb_pts{ind,1} = Pts{iGb};
                    
                    tb_iE_inGrain(ind,1) = iE_inGrain;
                    tb_iE_atBoundary(ind,1) = iE_atBoundary;
                    tboundary(uniqueBoundary == GBs(iGb)) = 1;
                else
                    disp('error');
                end
            end
        end
    end
    
    vN = {'iE','gb','g','ts','iE_g','iE_gb'};
    % TA = cell2table(cell(0,length(vN)));
    TA = array2table([tb_iE, tb_gbNum, tb_gNum, tb_tsNum, tb_iE_inGrain, tb_iE_atBoundary]);
    TA.Properties.VariableNames = vN;
    
    % plot, be careful
    tBoundaryCell{iE} = [];
    tBoundaryCell{iE} = tboundary;
end

iE = iE_current;    % change back


%%
timeStr = datestr(now,'yyyymmdd_HHMM');
try
    save(fullfile(saveDataPath, [sampleName,'_5_twin_result_with_gb_intersection_manual.mat']), 'TA', 'tb_gbNum', 'tb_iE', 'tb_gNum', 'tb_tsNum', 'tb_pts', 'tb_iE_inGrain', 'tb_iE_atBoundary',...
        'tBoundaryCell','tBoundary_accum','struCell','-append');
    % save([timeStr,'_twin_at_boundary_result_ws.mat'],'-v7.3');
catch
    save(fullfile(saveDataPath, [sampleName,'_5_twin_result_with_gb_intersection_manual.mat']), 'TA', 'tb_gbNum', 'tb_iE', 'tb_gNum', 'tb_tsNum', 'tb_pts', 'tb_iE_inGrain', 'tb_iE_atBoundary',...
        'tBoundaryCell','struCell','-v7.3');
end




