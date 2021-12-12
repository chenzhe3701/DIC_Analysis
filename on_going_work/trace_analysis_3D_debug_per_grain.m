% chenzhe, 2018-09-19
%
% After running 3D_(), need some modification.
% For correction, select a grain, and correct all strain levels.
%
% chenzhe, 2018-10-20.
% This is the debug version of analysis.

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
load('cityDistMap.mat');

% store all the clusterNumMap s, omit stop-0
clusterNumberMapCell = cell(1,length(STOP)-1);
struCell = cell(1,length(STOP)-1);
for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru','clusterNumMap','clusterNumMapCleaned');
    clusterNumberMapCell{iE} = clusterNumMapCleaned;
    
end
% load labeled results for debugging: 
load('twinMaps.mat');

%% plot something of interest
close all;
iE_select = 2;
try
    myplot(X,Y,twinMapCell{iE_select},boundaryTFB);caxis([18,24]);
catch
end
myplot(X,Y,clusterNumberMapCell{iE_select},boundaryTFB);

%% get ID from map, then
ids = find_ID_on_map(X,Y,ID,gcf,gca);   % 190

%% plot unit cell to check trace, one at a time
ind = find(gID==ids(1));
hcp_cell('euler',[gPhi1(ind),gPhi(ind),gPhi2(ind)], 'ss', 25:30, 'stress', [-1 0 0; 0 0 0; 0 0 0]);

%% (1) analyze
% get the local indices

stru = struCell{iE_start};
% for iS = 1 %1:length(stru)
% 181, 262, 1350, 1390, 697, 193, 442, 1016, 1532
iS = find(arrayfun(@(x) x.gID == ids(1),stru));  % for debugging. [for WE43, some grains: 378, 694, 1144] [697 interesting as there is a non-twin trace], 

% clean old data to relabel.
for iE = iE_start:iE_stop
    struCell{iE}(iS).cActiveSS = zeros(size(struCell{iE}(iS).cActiveSS));
end

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

% initialize for each grain (iS)
for iE = iE_start:iE_stop
    twinMapLocal{iE} = zeros(size(ID_local));
    sfMapLocal{iE} = zeros(size(ID_local));
end
twinMapCell_cluster = [];
sfMapCell_cluster = [];

% strain level at which the analysis was started -- it enters at every possible strain level, starting at iE_start 
iLoop_iE = iE_start;     
iLoop_iC = 1;
iLoop_iEC = 1;
cVolPctOld = 0;
cVolPctNotDecrease = 1;

%% This is the loop to run
close all;
% for iE_outer = iE_start:iE_stop
iE_entry = iLoop_iE;

% for each iE_outer, iC_outer, find the tracked iE_list, iC_list.

% for iC_outer = 1:length(struCell{iE_outer}(iS).cLabel)
iC_entry = iLoop_iC;

[iE_list, iC_list] = find_tracked_iE_iC_list(struCell, iS, iE_entry, iC_entry);

% Analyze all the linked iEs.  So, if iE_list(1)==iE_outer, it means it has not been analyzed before, then do [iE_list(ii),iC_list(ii)] pairs

iEC = iLoop_iEC;

% for iEC = 1:length(iE_list)

iE = iE_list(iEC);
iC = iC_list(iEC);
disp(['------------------------ [iE_outer, iC_outer, iE, iC] = [',num2str(iE_entry),', ',num2str(iC_entry),', ',num2str(iE),', ',num2str(iC),'] ------------------------']);

if iE_list(1) == iE_entry
    close all;
    clusterNumMapL = clusterNumberMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    clusterNumMapL(ID_local~=ID_current) = 0;  % First, clean-up those doesn't belong to this grain
    if debugTF
        myplot(clusterNumMapL);
    end
    
    ssAllowed = ones(ntwin,1);
    [twinMapCell_cluster, sfMapCell_cluster, struCell, haveActiveSS] = label_twin_trace(twinMapCell_cluster, sfMapCell_cluster, clusterNumberMapCell,x_local,y_local, indR_min,indR_max, indC_min,indC_max, ID_local,ID_current,...
        struCell,iS,iE,iC,iE_list,iC_list,iEC,iE_stop,traceND,traceSF,sampleMaterial,'twin',debugTF+1, 0.3,0.3,ssAllowed);
                    
%     twinMapLocal{iE} = twinMapLocal{iE} + fragments;
    
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
if iLoop_iC > length(struCell{iE_entry}(iS).cLabel)
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
disp(twinMapCell_cluster);
%% update twinMapLocal with tMapCell
for ii = 1:size(twinMapCell_cluster,1)
   for jj = 1:size(twinMapCell_cluster,2)
      if ~isempty(twinMapCell_cluster{ii,jj}) 
         twinMapLocal{ii} = twinMapLocal{ii} + twinMapCell_cluster{ii,jj}; 
         sfMapLocal{ii} = sfMapLocal{ii} + sfMapCell_cluster{ii,jj}; 
      end
   end
   if ~isempty(twinMapLocal{ii})
       myplot(twinMapLocal{ii}); title(['twinMap at strain: ',num2str(ii)]);
       myplot(sfMapLocal{ii}); title(['sfMap at strain: ',num2str(ii)])
   end
end
%% update. First clean old map, then add new map.
for iE = iE_start:iE_stop
    %     toClean = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    %     toClean(ID_local ~= ID_current) = 0;
    %     twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) - toClean + twinMapLocal{iE};
    %
    %     toClean = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    %     toClean(ID_local ~= ID_current) = 0;
    %     sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + sfMapLocal{iE};
    
    map_local = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);  % (1) Cut a squared map from big map
    twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;          % (2) Eliminate this sqaured region from the big map
    map_local(ID_local == ID_current) = 0;                              % (3) clean the grain area on the cut map
    map_local = map_local + twinMapLocal{iE};                           % (4) update the grain area on the cut map
    twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = twinMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;  % (5) Add the modified cut map to big map
    
    map_local = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;
    map_local(ID_local ~= ID_current) = 0;
    map_local = map_local + sfMapLocal{iE};
    sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = sfMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;
    
    map_local = cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
    cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = 0;
    map_local(ID_local ~= ID_current) = 0;
    map_local = map_local + cToGbDistMapLocal{iE};
    cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) = cToGbDistMapCell{iE}(indR_min:indR_max, indC_min:indC_max) + map_local;
    
end


%%

timeStr = datestr(now,'yyyymmdd_HHMM');
save([timeStr,'_twinMaps.mat'],'twinMapCell','-v7.3');


