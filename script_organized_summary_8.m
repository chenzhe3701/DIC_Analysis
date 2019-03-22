% script for analyzing twin-grain boundary intersection

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

% load previous twin_gb interaction result
load(fullfile(twinGbIntersectionPath, twinGbIntersectionFile));

%% (8) For each twin, look at its contacting neighbor, look at the strain components of twin in the neighbor grain 
% 

warning off;
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);
gamma = 0.1289; % twin shear for Mg

% For crystal coord: x = a_dir = [11-20], y = [1-100], z = [0001]
% Physical meaning of eij components: slip in i-direction, on plane with j-normal
% exy: single prismatic slip
% exz: single basal slip
% eyx: double prismatic slip
% eyz: double basal slip
% ezx, ezy: twinning or <c+a> slip
% [~body, single_prism, single_basal;]
% [double prism, ~body, double_basal;]
% [twin_<c+a>, twin_<c+a>, ~body]



% For iE = 2:5, for each grain, find accumulated active TS, and each TS's related neighbor
for iE = 4
    T = cell2table(cell(0,8));
    T.Properties.VariableNames={'ID','ts','sf','gb','activeTF','dp','db','dt'};
    for iS = 1:length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        
        ind = find(gID==ID_current);
        euler_1 = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        if (1==eulerAligned)
            g_1 = euler_to_transformation(euler_1,[0,0,0],[0,0,0]);
        else
            g_1 = euler_to_transformation(euler_1,[-90,180,0],[0,0,0]); % setting-2
        end
        
        all_neighbor = gNeighbors(ind,1:gNNeighbors(ind));   % [var] all neighbors, regardless of touching twin or not
        
        % get accumulative tGB from iE=2 to iE=current iE. Note that this 'unique' operation is not generally applicable
        tGb = [];   all_tGb = [];
        for ii = 1:6
            tGb{ii} = [];
        end
        for ie = 2:iE
            for ii=1:6
                tGb{ii} = unique([tGb{ii}, struCell{ie}(iS).tGb{ii}]);
                all_tGb = [all_tGb, tGb{ii}];
            end
        end
        all_tGb = unique(all_tGb);  % [var] find all grain boudnaryies of this matrix
        
        all_twin_neighbor = [];     % [var] find all twin-related neighbors
        for igb = 1:length(all_tGb)
            uniqueGB = all_tGb(igb);
            % find nb grain id
            gPair = [floor(uniqueGB/10000), mod(uniqueGB,10000)];
            ID_neighbor = gPair(~ismember(gPair, ID_current));
            all_twin_neighbor = [all_twin_neighbor, ID_neighbor];
        end
        
        for it = 1:6
            sf = struCell{iE}(iS).tSF(it);
            gbs_considered = tGb{it};    % the unique boundaries the twin intersect, may be empty
            twinned_TF = true;
            
            if isempty(gbs_considered)
                % reconstruct a 'gb_considered'
                gbs_considered = [10000*all_neighbor(all_neighbor>ID_current) + ID_current, 10000*ID_current + all_neighbor(all_neighbor<ID_current)];
                gbs_considered = [10000*all_twin_neighbor(all_twin_neighbor>ID_current) + ID_current, 10000*ID_current + all_twin_neighbor(all_twin_neighbor<ID_current)];
                twinned_TF = false;
            end
            
            % if twinned, considering the impinging grain boundary, if not, considering all /(all twin related) neighbors 
            for igb = 1:length(gbs_considered)
                uniqueGB = gbs_considered(igb);
                % find nb grain id
                gPair = [floor(uniqueGB/10000), mod(uniqueGB,10000)];
                ID_neighbor = gPair(~ismember(gPair, ID_current));
                
                ind = find(gID==ID_neighbor);
                euler_2 = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
                if (1==eulerAligned)
                    g_2 = euler_to_transformation(euler_2,[0,0,0],[0,0,0]);
                else
                    g_2 = euler_to_transformation(euler_2,[-90,180,0],[0,0,0]); % setting-2
                end
                g = g_2 * g_1'; % from grain_1 to grain_2
                
                iss = it + nss * 0;   % for Mg
                n = ss(1,:,iss);
                b = ss(2,:,iss);
                dispGrad = gamma * b' * n;   % displacement gradient, in crystal coordinate
                F = eye(3) + dispGrad;
                epsilon = (F'*F-eye(3))/2;              % in crystal coordinate
                
                % (1) expressed in neighbor
                dgm = g * dispGrad * g';   % disGrad expressed in grain_2 coordinate
                epsilonn = ((eye(3)+dgm)'*(eye(3)+dgm)-eye(3))/2;
                dp = abs(dgm(2))+abs(dgm(4));   % an invariant representing the shear on prismatic
                db = sqrt(dgm(7)^2+dgm(8)^2);   % an invariant representing the shear on basal
                dt = sqrt(dgm(3)^2+dgm(6)^2);   % an invariant representing the shear by twinning/<c+a> slip
                
                % --> insert codes here to debug
                
                T = [T; {ID_current, it+nss, sf, uniqueGB, twinned_TF, dp, db, dt}];
            end
            
        end
        
    end
end


warning on;

%%
ind = T.sf>0.3;
t = T(ind,:);
grpstats(t,'activeTF')


%%
% --> recorded codes for debug

%                 norm(epsilonn(:))
%                 abs(dgm(2))+abs(dgm(4))
% %                 dgm(2)^2+dgm(4)^2
% %                 abs(dgm(7))+abs(dgm(8))
%                 dgm(7)^2+dgm(8)^2
% %                 abs(dgm(3))+abs(dgm(6))
%                 dgm(3)^2+dgm(6)^2
%
%                 % (1.2) equivalent - 1
%                 disp('-------------------- 1');
%                 gg = euler_to_transformation([60,0,0],[0 0 0],[0 0 0]);
%                 dgm_1 = gg * dgm * gg'
%                 epsilon_1 = ((eye(3)+dgm_1)'*(eye(3)+dgm_1)-eye(3))/2
%
%                 norm(epsilon_1(:))
%                 abs(dgm_1(2))+abs(dgm_1(4))
% %                 dgm_1(2)^2+dgm_1(4)^2
% %                 abs(dgm_1(7))+abs(dgm_1(8))
%                 dgm_1(7)^2+dgm_1(8)^2
% %                 abs(dgm_1(3))+abs(dgm_1(6))
%                 dgm_1(3)^2+dgm_1(6)^2
%
%
%                 % (1.3) equivalent - 2
%                 disp('-------------------- 2');
%                 gg = euler_to_transformation([120,0,0],[0 0 0],[0 0 0]);
%                 dgm_2 = gg * dgm * gg'
%                 epsilon_2 = ((eye(3)+dgm_2)'*(eye(3)+dgm_2)-eye(3))/2
%
%                 norm(epsilon_2(:))
%                 abs(dgm_2(2))+abs(dgm_2(4))
% %                 dgm_2(2)^2+dgm_2(4)^2
% %                 abs(dgm_2(7))+abs(dgm_2(8))
%                 dgm_2(7)^2+dgm_2(8)^2
% %                 abs(dgm_2(3))+abs(dgm_2(6))
%                 dgm_2(3)^2+dgm_2(6)^2