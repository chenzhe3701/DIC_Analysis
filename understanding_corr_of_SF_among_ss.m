
% chenzhe, 2018-01-17
% randomly generate lots of euler angles, show the correlation between
% schmid factors of different slip systems.

[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Mg','twin');
ss = crystal_to_cart_ss(ssa,c_a);
stressTensor = [-1 0 0; 0 0 0; 0 0 0];
for i=1:1000
    euler = rand(1,3)*360;
    g = euler_to_transformation(euler,[0 0 0],[0 0 0]);
    for iss=1:24   % for Mg
        %         disp('---');
        N(iss,:) = ss(1,:,iss) * g;
        M(iss,:) = ss(2,:,iss) * g;
        SFs(i,iss) = N(iss,:) * stressTensor * M(iss,:)';
    end
end
figure; scatter(SFs(:,1),SFs(:,2)); title('x: SF ss# 1, y: SF ss# 2');
figure; scatter(SFs(:,1),SFs(:,1)+SFs(:,2)+SFs(:,3));   title('x: SF ss# 1, y: SF ss# 1+2+3');
figure; scatter(SFs(:,19),SFs(:,22));   title('x: SF ss# 19, y: SF ss# 22');