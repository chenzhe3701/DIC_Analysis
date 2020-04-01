

% chenzhe, 2018-05-10, randomly generate some euler angles.
% Show that a pair of twin systems always show correlated strain components
%
% chenzhe, 2018-05-28, add effective strain

clear;
clc;
close all; 

n = 1000;
[ss,c_a,ssa] = define_SS_cart('Mg','twin');

% initialize
exxPred = zeros(n,10);    % [gID, predicted_strain_of_twin_systems]
exxPred(:,1) = (1:n)';
exyPred = exxPred;
eyyPred = exxPred;
eeffPred = exxPred;
eeff3Pred = exxPred;
hWaitbar = waitbar(0,'finding twin region for grains ...');
for iS = 1:n
    
    euler = randi(360,1,3);  % or, can simulate
    
    g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
  
    gamma = 0.1289; % twin shear for Mg
    for iss=19:24   % for Mg
        %         disp('---');
        N(iss,:) = ss(1,:,iss) * g;
        M(iss,:) = ss(2,:,iss) * g;
        MN2{iss} = M(iss,:)'*N(iss,:);
        
        F3 = eye(3) + gamma*MN2{iss};
        epsilon3 = (F3'*F3-eye(3))/2;
        eeff3Pred(iS,iss-17) =  sqrt(2/3*(epsilon3(1).^2 + epsilon3(5).^2 + epsilon3(9).^2  + 2*epsilon3(2).^2 + 2*epsilon3(3).^2 + 2*epsilon3(6).^2 ));
        
        MN2{iss} = MN2{iss}(1:2,1:2);
        F = eye(2) + gamma*MN2{iss};
        epsilon = (F'*F-eye(2))/2;
        
        exxPred(iS,iss-17) = epsilon(1);    % [19:24] - 17 = [2:7]
        exyPred(iS,iss-17) = epsilon(2);
        eyyPred(iS,iss-17) = epsilon(4);
        eeffPred(iS,iss-17) = calculate_effective_strain(epsilon(1),epsilon(2),epsilon(4));

    end    

    waitbar(iS/n, hWaitbar);
end
close(hWaitbar);
%%
figure; plot(exxPred(:,2),exxPred(:,5),'.'); 
title('\epsilon_{xx}','fontweight','normal');
xlabel('Twin System # 1');
ylabel('Twin System # 4');
set(gca,'fontsize',18,'xlim',[-0.08 0.08],'ylim',[-0.08,0.08]);
axis square;
print('twin pair exx.tif','-dtiff');

figure; plot(exyPred(:,2),exyPred(:,5),'.'); 
title('\epsilon_{xy}','fontweight','normal');
xlabel('Twin System # 1');
ylabel('Twin System # 4');
set(gca,'fontsize',18,'xlim',[-0.08 0.08],'ylim',[-0.08,0.08]);
axis square;
print('twin pair exy.tif','-dtiff');

figure; plot(eyyPred(:,2),eyyPred(:,5),'.'); 
title('\epsilon_{yy}','fontweight','normal');
xlabel('Twin System # 1');
ylabel('Twin System # 4');
set(gca,'fontsize',18,'xlim',[-0.08 0.08],'ylim',[-0.08,0.08]);
axis square;
print('twin pair eyy.tif','-dtiff');

%% This demonstrate that, if 3D is considered, effective strain is always 0.0747
figure; histogram(eeffPred(:,2));
figure; histogram(eeffPred(:,3));
figure; histogram(eeffPred(:,5));
figure; histogram(eeff3Pred(:,2));
figure; histogram(eeff3Pred(:,3));
figure; histogram(eeff3Pred(:,5));


% figure; scatter(exxPred(:,3),exxPred(:,6)); title('exx ts-20 vs ts-23');
% figure; scatter(exxPred(:,4),exxPred(:,7)); title('exx ts-21 vs ts-24');
% figure; scatter(exyPred(:,3),exyPred(:,6)); title('exy ts-20 vs ts-23');
% figure; scatter(exyPred(:,4),exyPred(:,7)); title('exy ts-21 vs ts-24');
% figure; scatter(eyyPred(:,3),eyyPred(:,6)); title('exx ts-20 vs ts-23');
% figure; scatter(eyyPred(:,4),eyyPred(:,7)); title('exx ts-21 vs ts-24');
