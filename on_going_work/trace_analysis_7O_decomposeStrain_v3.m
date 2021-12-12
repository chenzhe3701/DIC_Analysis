% Chenzhe 2016-3-8
%
% For a specific slip system, write out the deformation gradient F.  Then
% calculate the predicted strain E = 1/2(F'*F-I).  F is proportional to the
% shear increment of a deformation system, but no need to be infinitesimal,
% because the form of strain is finite strain.
%
% Fit the strain point-wise for, e.g., grain(125)
%
% chenzhe, 2017-06-08.  Review codes.  Add fcc.
%
% note, required functions in chenFunctions
% ---------- level 1----------
% define_SS()
% crystal_to_cart_ss()
% euler_to_transformatoin()
%
% try to used method_5, directly use MNij and NMij

addChenFunction;

% traceStruct is for finding out traces
[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;

directory_DIC = uigetdir('','pick base DIC directory');
pathSave = uigetdir('','pick a directory to save the processed data');

STOP = {'001','002','003','004','005','006','007'};
STOP = {'0','1','2','3','4','5','6','7'};
iE_start = 2;   % elongation levels to analyze. 1-based. ------------------------------------ modify this -----------------
iE_stop = length(STOP);
% file name prefixes ---------------------------------------------------------------- could modify  ---------------
f1 = 'T5_#7_stop_';    %'20170430_ts5Al_02_e';
f2 = '_';

previousData = 'WE43_T6_C1_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data. --------------------
load(previousData);
% sampleName='Ti7Al_B2'; sampleMaterial='Ti'; stressTensor=[1 0 0; 0 0 0; 0 0 0]; % These are some of the prevously loaded settings.


TARGETIDS = [191]; % Mg: [89,129,135,191,201,210,327,401,422,553] from 3/4/5-7, Ti: [151,154,173,175,207,247], from stop 5-10

%% info is a 3D-matrix for the data interested at each stop(page)
IDs = arrayfun(@(x) x.ID, traceStruct); 
if(isempty(TARGETIDS))
    TARGETIDS=IDs; 
end

% stress = [0,200,400,500,550, 600,650,700,750,780, 790,800,810,810,810];
% stressTensor = [1 0 0; 0 0 0; 0 0 0];
[ssa,c_a] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);
nSS = size(ssa,3);

for iE = 1:iE_stop     % selected stops -------------------------------------------------------------
    
    strainFile = [directory_DIC,'\',f2,STOP{iE}]; disp(strainFile);    % note/change the prefix of DIC file %%%%%%%%%%%%%%%%%%%%
    loaded = load(strainFile,'exx','exy','eyy','sigma','exy_corrected');    % load strain data
    exx = loaded.exx;
    exy = loaded.exy;
    eyy = loaded.eyy;
    sigma = loaded.sigma;
    clear loaded;
    
    for ID_1 = TARGETIDS  % selected_ID -----------------------------------------------------------------------------------
        disp(['File: ',num2str(iE), ', Grain: ',num2str(ID_1)]);
        % (0) Extract local xyuv
        ind_gID = find(ID_1 == gID);    % an index of row
        
        ind_pool = ismember(ID, ID_1);
        indC_min = find(sum(ind_pool, 1), 1, 'first');
        indC_max = find(sum(ind_pool, 1), 1, 'last');
        indR_min = find(sum(ind_pool, 2), 1, 'first');
        indR_max = find(sum(ind_pool, 2), 1, 'last');
        
        nRow = indR_max - indR_min + 1;
        nColumn = indC_max - indC_min + 1;
        
        exx_current = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor
        exy_current = exy(indR_min:indR_max, indC_min:indC_max);
        eyy_current = eyy(indR_min:indR_max, indC_min:indC_max);
        sigma_current = sigma(indR_min:indR_max, indC_min:indC_max);
        
%         boundaryTF_current = boundaryTF{1}(indR_min:indR_max, indC_min:indC_max);
%         x_current = X(indR_min:indR_max, indC_min:indC_max);
%         y_current = Y(indR_min:indR_max, indC_min:indC_max);
        ID_map_current = ID(indR_min:indR_max, indC_min:indC_max);
        
        ind_not_grain = (ID_map_current~=ID_1);
        exx_grain = exx_current;
        exx_grain(ind_not_grain) = nan;  % 'e_grain' is strain of This grain. 'e_current' is strian of this region.
        exy_grain = exy_current;
        exy_grain(ind_not_grain) = nan;
        eyy_grain = eyy_current;
        eyy_grain(ind_not_grain) = nan;
%         
%         x_grain = x_current;
%         y_grain = y_current;
        sigma_grain = sigma_current;
        sigma_grain(ind_not_grain) = -1;
        
        ID_grain = ID_map_current;
        ID_grain(ind_not_grain) = nan;
        
        % myplot(exx_grain);caxis([-0.02 0.1]);title('\epsilon_x_x','fontsize',18);
        
        % (1) calculate mij of all slip systems for this grain
        ind_euler = find(gID==ID_1);
        euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
        g = euler_to_transformation(euler,[-90,180,0],[180,0,0]);     % ------------------------ your sample is flipped !!! 
        mij = [];   % tensor product of M and N, for 24 slip systems
        Mij = [];
        MNij = [];
        for ii = 1:nSS
            N(ii,:) = ss(1,:,ii) * g;
            M(ii,:) = ss(2,:,ii) * g;
            mij(:,:,ii) = (M(ii,:)'* N(ii,:) + N(ii,:)'*M(ii,:) )/2;
            Mij(:,:,ii) = N(ii,:)'*M(ii,:) * M(ii,:)'*N(ii,:)/2;
            MNij(:,:,ii) = M(ii,:)'* N(ii,:);
            NMij(:,:,ii) = N(ii,:)'* M(ii,:);
        end
        
        %     % (2) calculate Elastic strain of this grain
        %     [~,S] = elastic_constant_Ti();
        %     S = S/10^6; % Make it to MPa^(-1)
        %
        %     stressTensor_1 = stressTensor * stress(iFile);
        %     sElasticR = g*stressTensor_1*g';
        %     sElasticVector = sElasticR([1,5,9,6,3,2]);
        %     eElasticVector = S*sElasticVector';
        %     eElasticR = reshape(eElasticVector([1,6,5,6,2,4,5,4,3]),3,3).*[1, 1/2, 1/2; 1/2, 1, 1/2; 1/2, 1/2, 1];
        %     eElastic = g'*eElasticR*g;
        
        % (3) For each point on the [high resolution] strain map, find the best
        % fit dg using all active slip system information.
        
        ind = find(ID_1 == IDs);
        activeSSS = find(traceStruct(ind).ssActivated);       % active Slip SystemS
        activeSSS = [20,23];     % selected slip system -------------------------------------------------------------- manually
        
        % chenzhe, 2017-06-08.
        % For Al, always allow the 3 on the same plane to be active.
        if any(strcmpi(sampleMaterial, {'Al','Aluminum'}))
            expanded = [];
            if ~isempty(intersect(activeSSS,[1,2,3]))
                expanded = [expanded, 1,2,3];
            end
            if ~isempty(intersect(activeSSS,[4,5,6]))
                expanded = [expanded, 4,5,6];
            end
            if ~isempty(intersect(activeSSS,[7,8,9]))
                expanded = [expanded, 7,8,9];
            end
            if ~isempty(intersect(activeSSS,[10,11,12]))
                expanded = [expanded, 10,11,12];
            end
            activeSSS = expanded;
        end

        m11 = reshape(mij(1,1,:),nSS,1);
        m12 = reshape(mij(1,2,:),nSS,1);
        m22 = reshape(mij(2,2,:),nSS,1);
        M11 = reshape(Mij(1,1,:),nSS,1);
        M12 = reshape(Mij(1,2,:),nSS,1);
        M22 = reshape(Mij(2,2,:),nSS,1);
        
        dg_iTF = zeros(1,1,nSS);
        activeTF = zeros(1,1,nSS);
        penalty = zeros(1,1,nSS);
        
        options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');
        
        trimmer = zeros(3,3,nSS);
        for activeSS = activeSSS     
            dg_iTF(activeSS) = 1;    % initialized guess of dg !
            activeTF(activeSS) = 1;
            trimmer(:,:,activeSS) = [1 1 0; 1 1 0; 0 0 0];
        end
                
        [nR,nC] = size(exx_grain);
        dgCell = cell(size(exx_grain));
        costMat = zeros(size(exx_grain))*nan;
        flagMat = zeros(size(exx_grain))*nan;
        
        dgWeight = ones(size(exx_grain,1),size(exx_grain,2),nSS);
        if (iE > iE_start)&&exist([pathSave,'\slip_decompose_stop_',num2str(iE-1),'_grain_',num2str(ID_1)])
            loaded = load([pathSave,'\slip_decompose_stop_',num2str(iE-1),'_grain_',num2str(ID_1)],'dgCell');
            dgCellLoaded = loaded.dgCell;
            mFilter = fspecial('average',21);
            for activeSS = activeSSS
                if length(dgCellLoaded{activeSS})>1
                    weight = filter2(mFilter,dgCellLoaded{activeSS},'same');
                    weight = abs(weight);
                    baseLine = quantile(weight(:),0.2);
                    gap = quantile(weight(:),0.9) - baseLine;
                    weight = (weight - baseLine)/(0.5*gap);
                    weight = 1 + logsig(weight-1);
                    dgWeight(:,:,activeSS) = weight;
                end
            end
            clear loaded; clear dgCellLoaded;
        end
        maxWeight = max(dgWeight,[],3);
        for activeSS = activeSSS
            dgWeight(:,:,activeSS) =  maxWeight./dgWeight(:,:,activeSS);
        end
        
        % based on costFunction, determine best initial value to use.
        dg_pool = 4*(rand(1,1,nSS,20)-0.5);
        step = 40;  % do a rough check, real fitting use step = 5, here we use step = 40 and run 20 times, so add 5*5/40/40*20 = 30% more time 
        cost = zeros(1,20);
        for ii = 1:20
            dg_i = dg_pool(:,:,:,ii).*dg_iTF;
            parfor iR = 1:nR
                for iC = 1:nC
                    if (sigma_grain(iR,iC)~=-1)&&(rem(iR,step)==1)&&(rem(iC,step)==1)
                        exx_t = exx_grain(iR,iC);
                        exy_t = exy_grain(iR,iC);      % note the negative here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -------------------------------------- modify if needed
                        eyy_t = eyy_grain(iR,iC);
                        e_t = [exx_t, exy_t, 0; exy_t, eyy_t, 0; 0 0 0];  % note that exy should be reversed
                        penalty = zeros(1,1,nSS);
                        for activeSS = activeSSS
                            penalty(1,1,activeSS) = dgWeight(iR,iC,activeSS);
                        end
                        % e_t = [exx_t, exy_t, 0; exy_t, eyy_t, 0; 0 0 0];  % note that exy should be reversed
                        [dgCell{iR,iC},costMat(iR,iC),flagMat(iR,iC)] = fminunc(@(x) sum(sum((sum(1/2*(mtimes_page(eye(3)+NMij.*x, eye(3)+MNij.*x)-eye(3)).*trimmer,3)-e_t).^2)) ...
                                    + (sum(penalty.*sqrt(sum(sum( (1/2*(mtimes_page(eye(3)+NMij.*x, eye(3)+MNij.*x)-eye(3)).*trimmer).^2,1),2)),3)-sqrt(exx_t^2+eyy_t^2+2*exy_t^2)).^2 , dg_i, options);
                    else
                        dgCell{iR,iC} = zeros(1,1,nSS);
                    end
                end
            end
            cost(ii) = nanmean(costMat(:));
            disp(ii);
        end
        [~,indMin] = min(cost);
        dg_i = dg_pool(:,:,:,indMin).*dg_iTF;
        dgCell = cell(size(exx_grain));
        costMat = zeros(size(exx_grain))*nan;
        
        parfor_progress(nR);
        parfor iR = 1:nR
            for iC = 1:nC
                if (sigma_grain(iR,iC)~=-1)&&(rem(iR,5)==1)&&(rem(iC,5)==1)
                    exx_t = exx_grain(iR,iC);
                    exy_t = exy_grain(iR,iC);       % ---------------------------------------------------------------------- note here
                    eyy_t = eyy_grain(iR,iC);
                    e_t = [exx_t, exy_t, 0; exy_t, eyy_t, 0; 0 0 0];  % note that exy should be reversed
                    penalty = zeros(1,1,nSS);
                    for activeSS = activeSSS
                        penalty(1,1,activeSS) = dgWeight(iR,iC,activeSS);
                    end
                    % e_t = [exx_t, exy_t, 0; exy_t, eyy_t, 0; 0 0 0];  % note that exy should be reversed
                    [dgCell{iR,iC},costMat(iR,iC),flagMat(iR,iC)] = fminunc(@(x) sum(sum((sum(1/2*(mtimes_page(eye(3)+NMij.*x, eye(3)+MNij.*x)-eye(3)).*trimmer,3)-e_t).^2)) ...
                                    + (sum(penalty.*sqrt(sum(sum( (1/2*(mtimes_page(eye(3)+NMij.*x, eye(3)+MNij.*x)-eye(3)).*trimmer).^2,1),2)),3)-sqrt(exx_t^2+eyy_t^2+2*exy_t^2)).^2 , dg_i, options);
                    else
                    dgCell{iR,iC} = zeros(1,1,nSS);
                end
            end
            parfor_progress;
        end
        parfor_progress(0);
        
        dgMat = cell2mat(dgCell);
        clear dgCell;
        for ii=1:nSS
            if ismember(ii,activeSSS)
                dgCell{ii} = dgMat(:,:,ii);
                dgCell{ii}(dgCell{ii}==0) = nan;
                dgCell{ii}(ind_not_grain) = 0;
                dgCell{ii} = inpaint_nans(dgCell{ii},4);
                dgCell{ii}(ind_not_grain) = nan;
            else
                dgCell{ii} = 0;
            end
        end
        costMat = inpaint_nans(costMat,4);
        costMat(ind_not_grain) = nan;
        
        save([pathSave,'\slip_decompose_stop_',num2str(iE),'_grain_',num2str(ID_1)],'dgCell','costMat')
    
    end
    
end


