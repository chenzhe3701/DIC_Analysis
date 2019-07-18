% chenzhe, 2019-07-17
% It is likely useful to manually label twin-gb intersection.  We need some ground truth data anyway.    


clear;
addChenFunction;

% grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder to save the grain data'),'\'];
dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor','strainPauses');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path of the saved processed data, or WS, or etc.'),'\'];
saveDataPathInput = saveDataPath;
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
if ~strcmpi(saveDataPath,saveDataPathInput)
    disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
    return;
end

% Load from the pre-labeled results: twinMap, sfMap, struCell.  (cToGbDistMap is omitted, as will no longer be used in this code)
[confirmedLabelFile, confirmedLabelPath] = uigetfile('D:\p\m\DIC_Analysis\','select the results where twin identification was based on trace dir and strain');

[twinGbIntersectionFile, twinGbIntersectionPath] = uigetfile('D:\p\m\DIC_Analysis\20190222_1246_twin_at_boundary_result.mat','select the results for twin-grain boundary intersection');

try
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
catch
    load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted'],'X','Y','boundaryTF','boundaryTFB','uniqueBoundary','uniqueBoundaryList','ID','gID','gExx','gPhi1','gPhi','gPhi2','gNeighbors','gNNeighbors');
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





%%  manual label
continueTF = true;
iE=5;
iS = 162;
ID_current = 190;

struCell{iE}(iS).tGb = cell(1,ntwin);
struCell{iE}(iS).tGbPts = cell(1,ntwin);
struCell{iE}(iS).tGbNormal = cell(1,ntwin);
for ii = 1:ntwin
    struCell{iE}(iS).tGbPts{ii} = {[]};
    struCell{iE}(iS).tGbNormal{ii} = {[]};
end

while continueTF
    
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
    % mult_factor = ones(size(traceDir));
    % mult_factor(traceDir<0) = 1;
    % mult_factor(traceDir>=0) = -1;
    % traceND = traceDir + 90*mult_factor;    % convert traceDir to traceND
    
    
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
    boundaryTFB_local = boundaryTFB(indR_min:indR_max, indC_min:indC_max);
    trueTwinMapLocal = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    e_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
    
    [handleFig,aa,~] = myplot(X_local, Y_local, trueTwinMapLocal, boundaryTFB_local);
    
    % --> find existing lines in this region
    for iTwin = 1:6
        tGbNormals = struCell{iE}(iS).tGbNormal{iTwin};
        for iGb = 1:length(tGbNormals)
            if ~isempty(tGbNormals{iGb})
                for iLine = 1:length(tGbNormals{iGb})
                    drawline(aa,'Position',tGbNormals{iGb}{iLine});
                end
            end
        end
    end
    
    needToAdd = questdlg('add trace?','select answer','Yes','No','Cancel','No');
    switch needToAdd
        case 'Yes'
            nextTrace = true;
            while nextTrace
                % can use a button to start drawing a line
                handleDrawline = drawline(aa,'Color','r');
                posAdd = customWait(handleDrawline);
                
                answer = questdlg('add,redraw,or next grain?','select operation',...
                    'NextTrace','Redraw','AcceptAndNextGrain',...
                    'AcceptAndNextGrain');
                % find coordinate/maybe indices
                pt1 = posAdd(1,:);
                pt2 = posAdd(2,:);
                x = posAdd(:,1);
                y = posAdd(:,2);
                %             for ii = 1:size(x,1)
                %                 [~,subx] = min(abs(X_local(1,:)-x(ii)));
                %                 [~,suby] = min(abs(Y_local(:,1)-y(ii)));
                %                 ids(ii) = ID_local(suby,subx);
                %                 indr(ii) =  suby;
                %                 indc(ii) =  subx;
                %             end
                
                % find intersection between the line drawn and uniqueGrainBoundary.
                [pt,inds,indR,indC] = grids_covered_by_line(X_local,Y_local,pt1,pt2);
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
                    [pt,inds,indR,indC] = grids_covered_by_line(X_local,Y_local,pt1,pt2);
                    gbNum = uniqueBoundary_local(inds);
                    gbNum = gbNum(gbNum>0);
                end
                gbNum = gbNum(1);
                
                ind_of_ind = find(uniqueBoundary_local(inds)==gbNum, 1, 'first');
                ind = inds(ind_of_ind);
                xcoord = X_local(ind);
                ycoord = Y_local(ind);
                
                switch answer
                    case 'NextTrace'
                        % look at direction of line drawn, determine which twin system it is
                        direction = atand((y(2)-y(1))/(x(2)-x(1)));
                        [~,iTwin] = min(abs(direction-traceDir));
                        
                        ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                        if isempty(ind)
                            % append
                            struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                            iGb = length(struCell{iE}(iS).tGb{iTwin});
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                            struCell{iE}(iS).tGbNormal{iTwin}{iGb} = {[pt1;pt2]};
                        else
                            % add
                            iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} = [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                            struCell{iE}(iS).tGbNormal{iTwin}{iGb} = [struCell{iE}(iS).tGbNormal{iTwin}{iGb}; {[pt1;pt2]}];
                        end
                    case 'Redraw'
                        delete(handleDrawline);
                    case 'AcceptAndNextGrain'
                        % look at direction of line drawn, determine which twin system it is
                        direction = atand((y(2)-y(1))/(x(2)-x(1)));
                        [~,iTwin] = min(abs(direction-traceDir));
                        
                        ind = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                        if isempty(ind)
                            % append
                            struCell{iE}(iS).tGb{iTwin} = [struCell{iE}(iS).tGb{iTwin}, gbNum];   % append the gbNum of the gb touched by this twin
                            iGb = length(struCell{iE}(iS).tGb{iTwin});
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} =  [xcoord, ycoord];  % assign point coord (1x2 vector) to the cell value
                            struCell{iE}(iS).tGbNormal{iTwin}{iGb} = {[pt1;pt2]};
                        else
                            % add
                            iGb = find(struCell{iE}(iS).tGb{iTwin} == gbNum);
                            struCell{iE}(iS).tGbPts{iTwin}{iGb} = [struCell{iE}(iS).tGbPts{iTwin}{iGb}; [xcoord, ycoord]];    % append the point coord (1x2 vector) as new rows
                            struCell{iE}(iS).tGbNormal{iTwin}{iGb} = [struCell{iE}(iS).tGbNormal{iTwin}{iGb}; {[pt1;pt2]}];
                        end
                        close(handleFig);
                        nextTrace = false;
                end
            end
            
        case 'No'
            close(handleFig);
        case 'Cancle'
            close(handleFig);
            continueTF = false;
    end
    
    try
        save('struct_for_recover','struCell','-append');
    catch
        save('struct_for_recover','struCell');
    end
    
    disp(['iE=',num2str(iE),', iS=',num2str(iS)]);
end









