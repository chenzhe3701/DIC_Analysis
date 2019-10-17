% script for analyzing twin-grain boundary intersection
%
% chenzhe, 2019-06-26 note
% Similar to a ref paper, look at the displacement tensor components of
% twins in the neighbor grain's coordinate.
% And do some summary.

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
plotTF = 0;
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
for iE = 5
    T = cell2table(cell(0,16));
    T.Properties.VariableNames={'ID','ts','sf','tAF','gb',  'twinned_TF','twinned_nb_TF','twin_twin_gb_TF','db','dp',  'dt','basal_SF','prism_SF','twin_SF','miso', 'mPrime'};
    iS = 1;
    %%
    while iS<length(struCell{iE})
        ID_current = struCell{iE}(iS).gID;
        
        ind = find(gID==ID_current);
        euler_1 = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
        if (1==eulerAligned)
            g_1 = euler_to_transformation(euler_1,[0,0,0],[0,0,0]);
        else
            g_1 = euler_to_transformation(euler_1,[-90,180,0],[0,0,0]); % setting-2
        end
        
        % [1] [var] all neighbors, regardless of touching twin or not
        all_neighbor = gNeighbors(ind,1:gNNeighbors(ind));
        all_neighbor_copy = all_neighbor;  % copy for temporary use
        
        % [2.1] get accumulative tGB from iE=2 to iE=current iE. Note that this 'unique' operation is not generally applicable
        tGb = [];
        all_tGb = [];
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
        
        %         % [2.2] [var] then find all twin-related neighbors
        %         all_twin_neighbor = [];
        %         for igb = 1:length(all_tGb)
        %             uniqueGB = all_tGb(igb);
        %             % find nb grain id
        %             gPair = [floor(uniqueGB/10000), mod(uniqueGB,10000)];
        %             ID_neighbor = gPair(~ismember(gPair, ID_current));
        %             all_twin_neighbor = [all_twin_neighbor, ID_neighbor];
        %         end
        
        if ~isempty(all_tGb)
            activeTS = find(cellfun(@(x) ~isempty(x),tGb));
        end
                    
        if plotTF==1
            % ////// find out local data for local map, ID_local, e_local, uniqueGB_local, ... for potential plot
            ind_local = ismember(ID, [ID_current;all_neighbor(:)]); %ismember(ID, [ID_current,ID_neighbor]);
            % Make it one data point wider on each side
            indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
            indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
            indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
            indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);
            
            ID_local = ID(indR_min:indR_max, indC_min:indC_max);
            x_local = X(indR_min:indR_max, indC_min:indC_max);
            y_local = Y(indR_min:indR_max, indC_min:indC_max);
            uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
            % map_local = trueTwinMapCell{iE}(indR_min:indR_max, indC_min:indC_max);
            map_local = strainFile{iE}.exx(indR_min:indR_max, indC_min:indC_max);
            
            % if this grain has twin, and we want to plot
            h_f = myplot(x_local, y_local, map_local, uniqueBoundary_local);
            hold on;
            label_map_with_text(x_local, y_local, ID_local, h_f, 'target_ID',ID_current, 'color','r', 'text',['TS=', num2str(activeTS)]);

        end
        
        % [3] for each twin, analyze the potentially related neighbor/ or just analyze the twin related neighbor
        for iTwin = 1:6
            sf = struCell{iE}(iS).tSF(iTwin);
            gbs_considered = tGb{iTwin};    % the unique boundaries the twin intersect, may be empty
            twinned_TF = true;
            
            if isempty(gbs_considered)
                % reconstruct a 'gb_considered'
                gbs_considered = [10000*all_neighbor(all_neighbor>ID_current) + ID_current, 10000*ID_current + all_neighbor(all_neighbor<ID_current)];
                % gbs_considered = [10000*all_twin_neighbor(all_twin_neighbor>ID_current) + ID_current, 10000*ID_current + all_twin_neighbor(all_twin_neighbor<ID_current)];
                twinned_TF = false;
            end
            
            % The assumption is that, we always find a twin intersecting a grain boundary.  But here we can provide a double check
            if ~isempty(all_tGb) && (ismember(iTwin,activeTS) && (twinned_TF==false))
                disp('This does not match the assumption that a twin always have at least one intersecting grain boundary');
            end
            
            tAF = struCell{iE}(iS).tVol(iTwin)/struCell{iE}(iS).gVol;  % twin area fraction of grain size
            
            % [3] if twinned, considering the impinging grain boundary, if not, considering all /(all twin related) neighbors
            for igb = 1:length(gbs_considered)
                uniqueGB = gbs_considered(igb);
                % find nb grain id
                gPair = [floor(uniqueGB/10000), mod(uniqueGB,10000)];
                ID_neighbor = gPair(~ismember(gPair, ID_current));  % [var] ID of the grain on the other side of this twin-boundary
                
                % find nb euler angle, rotation matrix
                ind = find(gID==ID_neighbor);
                euler_2 = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
                if (1==eulerAligned)
                    g_2 = euler_to_transformation(euler_2,[0,0,0],[0,0,0]);
                else
                    g_2 = euler_to_transformation(euler_2,[-90,180,0],[0,0,0]); % setting-2
                end
                g = g_2 * g_1'; % from grain_1 to grain_2
                
                % calculate displacement gradient, ...
                iss = iTwin + nss;   % for Mg
                n = ss(1,:,iss);
                b = ss(2,:,iss);
                dispGrad = gamma * b' * n;   % displacement gradient, in crystal coordinate
                F = eye(3) + dispGrad;
                epsilon = (F'*F-eye(3))/2;              % in crystal coordinate
                
                % Displacement gradient expressed in neighbor, 'displacement gradient matrix' for 'dgm'
                dgm = g * dispGrad * g';   % disGrad expressed in grain_2 coordinate
                epsilonn = ((eye(3)+dgm)'*(eye(3)+dgm)-eye(3))/2;
                dp = abs(dgm(2))+abs(dgm(4));   % an invariant representing the shear on prismatic
                db = sqrt(dgm(7)^2+dgm(8)^2);   % an invariant representing the shear on basal
                dt = sqrt(dgm(3)^2+dgm(6)^2);   % an invariant representing the shear by twinning/<c+a> slip
                % If consider symmetry here, just checking
                %                 [~,gsym] = hcp_symmetry;
                %                 for isym = 1:size(gsym,3)
                %                     dgm_sym = gsym(:,:,isym) * dgm * gsym(:,:,isym)'
                %                     dp = abs(dgm(2))+abs(dgm(4))
                %                     db = sqrt(dgm(7)^2+dgm(8)^2)
                %                     dt = sqrt(dgm(3)^2+dgm(6)^2)
                %                 end
                
                % For this neighbor, there is a matrix describing how the eij of 6 twin systesm can be expressed in the 24 systems (5 modes) in this neighbor.
                exz_iTwin_jMode = calculate_exz(euler_1, euler_2,sampleMaterial);
                ez = exz_iTwin_jMode./max(exz_iTwin_jMode,[],1);    % because it is normalized to 1, it might be reasonable to use the value (rather than rank) to represent the easiness. 
                
                % ////// find how what is happening actually in neighbor? What is the active system in this neighbor?
                iS_nb = find(arrayfun(@(x) x.gID == ID_neighbor, struCell{iE}));
                twinned_nb_TF = false;
                twin_twin_gb_TF = false;
                
                % ////// Calculate Basal/Twin schmid factor in this neighbor grain
                [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler_2, [0 0 0], [0 0 0], [-1 0 0; 0 0 0; 0 0 0], sampleMaterial, 'twin');
                basal_SF = max(abs_schmid_factor(1:3,2));
                prism_SF = max(abs_schmid_factor(4:6,2));
                [twin_SF, thatTS] = max(abs_schmid_factor(19:24,2));
                
                % ////// Need to calculate something else: misorentation, mPrime, etc ...   
                [miso, ~] = calculate_misorientation_euler_d(euler_1, euler_2, 'HCP');
                [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler_1, euler_2, [-1 0 0; 0 0 0; 0 0 0], [0 1 0], sampleMaterial, 'twin');
                mPrime = mPrimeMatrixAbs(19:24,1:3);    % of interest
                mPrime = max(mPrime(iTwin,:)); % max of this twin w.r.t a basal slip   
                
                if ~isempty(iS_nb)
                    % (1) If twinned?
                    activeTS_nb = sum(struCell{iE}(iS_nb).cTrueTwin,1);
                    if any(activeTS_nb)
                        twinned_nb_TF = true;
                    end
                    % (2) if the twin intersect this boundary?
                    all_tGb_nb = [];
                    for ie = 2:iE
                        for ii=1:6
                            all_tGb_nb = [all_tGb_nb, struCell{ie}(iS_nb).tGb{ii}];
                        end
                    end
                    all_tGb_nb = unique(all_tGb_nb);  % [var] find all grain boudnaryies of the neighbor grain
                    if ismember(uniqueGB,all_tGb_nb)
                        % twins in neighbor intersect this grain boundary
                        twin_twin_gb_TF = true;
                    end
                end
                
                if twinned_TF || true
                    T = [T; {ID_current, iTwin+nss, sf, tAF, uniqueGB, twinned_TF, twinned_nb_TF, twin_twin_gb_TF, db, dp, dt, basal_SF, prism_SF, twin_SF, miso, mPrime}];
                end
                
                if ~isempty(all_tGb) && twinned_TF && plotTF
                    % we can plot the local data, to examine
                    textString = ['pr=',num2str(dp,3), newline,'ba=',num2str(db,3), newline,'tw=', num2str(dt,3)];
                    textString = ['ba =',num2str(ez(iTwin,1)), newline,'pr =',num2str(ez(iTwin,2)), newline,'tw =', num2str(ez(iTwin,5))];  % new code after 24 symmetry
                    label_map_with_text(x_local, y_local, ID_local, h_f, 'target_ID', ID_neighbor, 'color', 'k', 'text', textString);
                    all_neighbor_copy(all_neighbor_copy==ID_neighbor) = [];
                end
                
            end
            
        end
        
        if ~isempty(all_tGb) && plotTF
            label_map_with_text(x_local, y_local, ID_local, h_f, 'target_ID', all_neighbor_copy(:), 'color', 'y', 'text', []);
        end
        
        % This is for manual debug purpose. Use while loop.  If have an example to plot, pause the program.
        iS = iS+1
        if ~isempty(all_tGb) && plotTF
            break;
        end
        
        
    end
end


warning on;

%% Compare [dp/db/dt] between [twinned] vs. [not-twinned], using histogram  
close all;
ind = T.sf>0.15;
t = T(ind,:);
grpstats(t,'twinned_TF')

ind = t.twinned_TF == 1;
t1 = t(ind,:);
t2 = t(~ind,:);
figure; hold on;
histogram(t1.dt); 
histogram(t2.dt);
legend('twinned','not twinned');
title('twin and <c+a>')

%%
ind = t.twinned_TF==true;
figure; hold on;
histogram(t.db(ind));
histogram(t.db(~ind));


%% [analyze] For the twinned grains, Twin_Area_Fraction vs DispGrad_dp/db/dt  
close all;
ind = T.twinned_TF==1;
t = T(ind,:);
figure;
plot(t.dt, t.tAF, '.');
xlabel('dt');
ylabel('twin area fraction');





