

% chenzhe, 2018-01-17
% look at for all the grains, how predicted twin systems strain conponents are distributed.
% based on trace_analysis_2, need to preload some strain and grain data.

close all;
% initialize
exxPred = zeros(size(gIDwithTrace,1),10);    % [gID, predicted_strain_of_twin_systems]
exxPred(:,1) = gIDwithTrace;
exyPred = exxPred;
eyyPred = exxPred;

hWaitbar = waitbar(0,'finding twin region for grains ...');
for iS =1:length(gIDwithTrace)
    close all;
    % select the target grain
    ID_current = gIDwithTrace(iS); 
    
    % ================ method-1, from theoretical twin shear --> to predict strain components as cluster center ===================
    ind_euler = find(gID==ID_current);
    euler = [gPhi1(ind_euler),gPhi(ind_euler),gPhi2(ind_euler)];
    if (1==eulerAligned)
        g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
    else
        g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
    end    
    gamma = 0.1289; % twin shear for Mg
    for iss=19:24   % for Mg
        %         disp('---');
        N(iss,:) = ss(1,:,iss) * g;
        M(iss,:) = ss(2,:,iss) * g;
        MN2{iss} = M(iss,:)'*N(iss,:);
        MN2{iss} = MN2{iss}(1:2,1:2);
        F = eye(2) + gamma*MN2{iss};
        epsilon = (F'*F-eye(2))/2;
        
        exxPred(iS,iss-17) = epsilon(1);    % [19:24] - 17 = [2:7]
        exyPred(iS,iss-17) = epsilon(2);
        eyyPred(iS,iss-17) = epsilon(4);

    end    

    waitbar(iS/length(gIDwithTrace), hWaitbar);
end
close(hWaitbar);

figure; scatter(exxPred(:,2),exxPred(:,5)); title('exx ts-19 vs ts-22');
figure; scatter(exxPred(:,3),exxPred(:,6)); title('exx ts-20 vs ts-23');
figure; scatter(exxPred(:,4),exxPred(:,7)); title('exx ts-21 vs ts-24');

figure; scatter(exyPred(:,2),exyPred(:,5)); title('exy ts-19 vs ts-22');
figure; scatter(exyPred(:,3),exyPred(:,6)); title('exy ts-20 vs ts-23');
figure; scatter(exyPred(:,4),exyPred(:,7)); title('exy ts-21 vs ts-24');

figure; scatter(eyyPred(:,2),eyyPred(:,5)); title('exx ts-19 vs ts-22');
figure; scatter(eyyPred(:,3),eyyPred(:,6)); title('exx ts-20 vs ts-23');
figure; scatter(eyyPred(:,4),eyyPred(:,7)); title('exx ts-21 vs ts-24');
