

% useful scripts:

%% (1) get unique boundary
[~, boundaryID, neighborID, ~, ~] = find_one_boundary_from_ID_matrix(ID);
uniqueBoundary = max(boundaryID,neighborID)*10000 + min(boundaryID,neighborID);
uniqueBoundaryList = unique(uniqueBoundary(:)); 
uniqueBoundaryList(uniqueBoundaryList==0) = [];

%% (2) select grain area

ID_current = 1026;
ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
% Make it one data point wider on each side
indC_min = max(1, find(sum(ind_local, 1), 1, 'first')-1);
indC_max = min(size(ID,2), find(sum(ind_local, 1), 1, 'last')+1);
indR_min = max(1, find(sum(ind_local, 2), 1, 'first')-1);
indR_max = min(size(ID,1), find(sum(ind_local, 2), 1, 'last')+1);

ID_local = ID(indR_min:indR_max, indC_min:indC_max);
uniqueBoundary_local = uniqueBoundary(indR_min:indR_max, indC_min:indC_max);
map_local = trueTwinMapCell{3}(indR_min:indR_max, indC_min:indC_max);
map_local = twinMapCell{3}(indR_min:indR_max, indC_min:indC_max);
%% get IS from grain ID
iS = find(arrayfun(@(x) x.gID == 164,struCell{iE}));

%% calculate traceSF and traceND
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

%% check if there are twins that are active in iE, but not active in iE-1
clc
for iE = 3:5
   for iS = 1:length(struCell{iE})
       activeTS = logical(sum(struCell{iE}(iS).cTrueTwin,1));
       activeTS_pre = logical(sum(struCell{iE-1}(iS).cTrueTwin,1));
       if any(activeTS-activeTS_pre<0)
           disp([iE,iS,struCell{iE}(iS).gID]);
           disp(activeTS_pre);
           disp(activeTS);
       end
   end    
end

%% check if struCell is modified
clc
for iE = 2:5
    for iS = 1:length(struCell{iE})
        if any(struCellA{iE}(iS).cTrueTwin(:) - struCell{iE}(iS).cTrueTwin(:))
           disp([iE,iS,struCell{iE}(iS).gID]);
           disp(struCellA{iE}(iS).cTrueTwin);
           disp(struCell{iE}(iS).cTrueTwin);
        end        
    end    
end




%% scripts for predicted strain
[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);

gamma = 0.1289; % twin shear for Mg
cPred = nan*zeros(nss,5);   % [iss, SF, exx, exy, eyy]
for iss = (nss+1):(nss+ntwin)   % for Mg
    %         disp('---');
    N(iss,:) = ss(1,:,iss) * g;
    M(iss,:) = ss(2,:,iss) * g;
    MN2{iss} = M(iss,:)'*N(iss,:);
    MN2{iss} = MN2{iss}(1:2,1:2);
    %         F3 = eye(3) + gamma*M(iss,:)'*N(iss,:);
    %         F = F3(1:2,1:2);
    F = eye(2) + gamma*MN2{iss};
    epsilon = (F'*F-eye(2))/2;
    %         disp((F3'*F3-eye(3))/2);
    %         disp(epsilon);
    cPred(iss,1) = iss;                                     % ss number
    cPred(iss,2) = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
    cPred(iss,3:5) = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.  Note that 'conjugated' twin system, i.e., 19 and 22, almost always show similar components!!!
end












