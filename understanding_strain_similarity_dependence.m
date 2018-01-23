
% Idea:
% The twin strain should be 0.1289 for Mg.  However, when identifying based
% measured strain components, a range of values should be accepted.
% So this code is intended to estimate the effect of variation in the 'twin
% shear' on the measured strain components.
%
% algorithm:
% (1) Randomly generate a euler angle, calculate the SF, and strain components.
% (2) Randomly alter the shear, calculate the corresponding strain components.
% (3) Calculate the distance between the unaltered and altered strain components.
% (4) Redistribute the difference between altered and unaltered strain
% components, and calculate a fitted shear, and the cost.
% Repeat for many euler angles, and many delta_gamma values.

[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Mg','twin');
ss = crystal_to_cart_ss(ssa,c_a);
stressTensor = [-1 0 0; 0 0 0; 0 0 0];

gamma = 0.1289;
delta_gamma = 0.1;
options = optimoptions(@fminunc,'display','off','algorithm','quasi-newton');
data = [];
for delta_gamma = -0.05:0.001:0.05
    for ii = 1:100
        euler = rand(1,3)*360;
        g = euler_to_transformation(euler,[0 0 0],[0 0 0]);
        
        for iss = nss+1:nss+ntwin
            N(iss,:) = ss(1,:,iss) * g;
            M(iss,:) = ss(2,:,iss) * g;
            MN2{iss} = M(iss,:)'*N(iss,:);
            MN2{iss} = MN2{iss}(1:2,1:2);
            
            F = eye(2) + gamma * MN2{iss};
            epsilon = (F'*F-eye(2))/2;
            SF = N(iss,:) * stressTensor * M(iss,:)';     % Schmid factor
            eij = [epsilon(1), epsilon(2), epsilon(4)];  % strain exx, exy, eyy.
            
            
            F_deviated = eye(2) + (gamma+delta_gamma) * MN2{iss};
            epsilon_deviated = (F_deviated'*F_deviated-eye(2))/2;
            eij_deviated = [epsilon_deviated(1), epsilon_deviated(2), epsilon_deviated(4)];  % strain exx, exy, eyy.
            
            d1 = pdist2(eij,eij_deviated,'cityblock');
            d2 = pdist2(eij,eij_deviated,'euclidean');
            d3 = pdist2(eij,eij_deviated,'squaredeuclidean');
            
            eij_shuffle = eij + (eij_deviated - eij) .* (rand(size(eij))-0.5)*2;
            
            [shear_fitted,cost]=fminunc(@(x) sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
                [eij_shuffle(1),eij_shuffle(2);eij_shuffle(2),eij_shuffle(3)]).^2)), 0.1298, options);
            
            %             % This is too much variation ...
            %             randomVec = rand(size(eij))-0.5;
            %             randomVec = randomVec./norm(randomVec,1);
            %             eij_shuffle_2 = eij + norm(eij_deviated - eij,1) .* (rand(size(eij))-0.5)*2;
            
            [shear_fitted_2,cost_2]=fminunc(@(x) sum(sum((((eye(2)+x*MN2{iss})'*(eye(2)+x*MN2{iss})-eye(2))/2-...
                [eij_shuffle_2(1),eij_shuffle_2(2);eij_shuffle_2(2),eij_shuffle_2(3)]).^2)), 0.1298, options);
            
            data = [data; delta_gamma, euler, eij, eij_deviated, SF, d1, d2, d3, shear_fitted, cost];
        end
    end
    disp(delta_gamma);
    % (1) delta_gamma, (2-4) euler, (5-7) eij, (8-10) eij_dev,
    % (11) SF, (12) d_city, (13) d_euc, (14) d_sqeuc, (15) gamma_fitted,
    % (16) cost,  % not implemented: (17) gamma_fitted_2, (18) cost_2
end
save('ShearTuneUpData.mat','data');
%% plot (0) 3 types of distances
figure;
plot3(data(:,12),data(:,13),data(:,14),'.');
xlabel('d-city'); ylabel('d-euclidean'); zlabel('d-sqeuclidean');
%% plot (1) delta_gamma vs d_euc, for different SFs
% for high SF ones, the d_euclidean is actually expected to be larger
ind_1 = data(:,11)>0.35;
ind_2 = ~ind_1;
figure;
subplot(1,2,1);
plot(data(ind_1,1),data(ind_1,13),'.');
title('SF > 0.3'); xlabel('delta-gamma'); ylabel('d-euclidean')
subplot(1,2,2);
plot(data(ind_2,1),data(ind_2,13),'.');
title('SF < 0.3'); xlabel('delta-gamma'); ylabel('d-euclidean')

%% plot (2) delta_gamma vs gamma_fitted vs cost, for different SFs
ind_1 = data(:,11)>0.35;
ind_2 = ~ind_1;
figure;
subplot(1,2,1);
plot3(data(ind_1,1),data(ind_1,15),data(ind_1,16),'.');
title('SF > 0.3'); xlabel('delta-gamma'); ylabel('fitted-gamma'); zlabel('cost')
subplot(1,2,2);
plot3(data(ind_2,1),data(ind_2,15),data(ind_2,16),'.');
title('SF < 0.3'); xlabel('delta-gamma'); ylabel('fitted-gamma'); zlabel('cost')

% use sqrt(cost)
figure;
subplot(1,2,1);
plot3(data(ind_1,1),data(ind_1,15),sqrt(data(ind_1,16)),'.');
title('SF > 0.3'); xlabel('delta-gamma'); ylabel('fitted-gamma'); zlabel('sqrt-cost')
subplot(1,2,2);
plot3(data(ind_2,1),data(ind_2,15),sqrt(data(ind_2,16)),'.');
title('SF < 0.3'); xlabel('delta-gamma'); ylabel('fitted-gamma'); zlabel('sqrt-cost')


%% Basically, the 'upper bound' should be more important than the'distribution'
% for a given tolerance in delta_gamma, use









