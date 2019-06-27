
% There needs to be an 'incoming' and 'outgoing' burgers vector.  Because
% the slip direction is bidirectional, the residual burgers vector will
% depend on the direction of the two burgers vectors.  Therefore, the
% direction of grain boundary should be described.  In fact the boundary
% should be specified everywhere, but here we just simplify it to be a
% straight line.
%
% Zhe Chen, 2016-1-18 
%
% function [schmidFactorG1, schmidFactorG2, mPrimeMatrixRaw,
% resBurgersMatrixRaw, mPrimeMatrixAbs, resBurgersMatrixAbs] =
% calculate_G1G2_mPrime_resB(phi1G1, phiG1, phi2G1, phi1G2, phiG2, phi2G2,
% phi1Sys, phiSys, phi2Sys, phiError, stressTensor, gbNormal)
%
% rewrite, 2019-04-02

function [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = ...
    calculate_mPrime_and_resB(euler_1, euler_2, stressTensor, gbNormal, material, twinTF)

m1 = angle2dcm(euler_1(1)/180*pi, euler_1(2)/180*pi, euler_1(3)/180*pi, 'zxz');
m2 = angle2dcm(euler_2(1)/180*pi, euler_2(2)/180*pi, euler_2(3)/180*pi, 'zxz');

% phiSys = -90, 180, 0
% phiError = 0;
% stressTensor = [1 0 0; 0 0 0; 0 0 0];

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(material, twinTF);
nSS = size(ssa,3);  %24
ss = crystal_to_cart_ss(ssa,c_a);

% Normalized slip direction and slip plane normal is used to calculate
% schmid factor.  But to calculate resBurgersVector, we need to consider
% that <c+a> slip is longer than <a> slip.
switch material
    case {'Mg','mg'}
        burgersWeight = ones(1,nSS);
        burgersWeight(13:18) = burgersWeight(13:18) * sqrt(1+c_a^2);
        burgersWeight(19:24) = burgersWeight(19:24) * sqrt(3+c_a^2);
end

for iSS = 1:nss
    slipPlaneG1(iSS,:) = ss(1,:,iSS) * m1;         % for grain 1, [slip plane normal] expressed in Global system
    slipDirectionG1(iSS,:) = ss(2,:,iSS) * m1;     % for grain 1, [slip direction] expressed in Global system
    slipDirectionG1(iSS,:) = slipDirectionG1(iSS,:) * sign(dot(slipDirectionG1(iSS,:), gbNormal));  % This makes the incoming and outgoing burgers vector pointing to the same side of the GB.
    slipDirectionG1W(iSS,:) = slipDirectionG1(iSS,:)*burgersWeight(iSS);
    schmidFactorG1(iSS,1) = abs(slipPlaneG1(iSS,:) * stressTensor * slipDirectionG1(iSS,:)');  % Schimid factor for slip system j
end
for iSS = nss+1 : nss+ntwin
    slipPlaneG1(iSS,:) = ss(1,:,iSS) * m1;         % for grain 1, [slip plane normal] expressed in Global system
    slipDirectionG1(iSS,:) = ss(2,:,iSS) * m1;     % for grain 1, [slip direction] expressed in Global system
    slipDirectionG1W(iSS,:) = slipDirectionG1(iSS,:)*burgersWeight(iSS);
    schmidFactorG1(iSS,1) =slipPlaneG1(iSS,:) * stressTensor * slipDirectionG1(iSS,:)';  % Schimid factor for TWIN system j
end



for iSS = 1:nss
    slipPlaneG2(iSS,:) = ss(1,:,iSS) * m2;         % for grain 1, [slip plane normal] expressed in Global system
    slipDirectionG2(iSS,:) = ss(2,:,iSS) * m2;     % for grain 1, [slip direction] expressed in Global system
    slipDirectionG2(iSS,:) = slipDirectionG2(iSS,:) * sign(dot(slipDirectionG2(iSS,:), gbNormal));      % This makes the incoming and outgoing burgers vector pointing to the same side of the GB.
    slipDirectionG2W(iSS,:) = slipDirectionG2(iSS,:)*burgersWeight(iSS);
    schmidFactorG2(iSS,1) = abs(slipPlaneG2(iSS,:) * stressTensor * slipDirectionG2(iSS,:)');  % Schimid factor for slip system j
end
for iSS = nss+1 : nss+ntwin
    slipPlaneG2(iSS,:) = ss(1,:,iSS) * m2;         % for grain 1, [slip plane normal] expressed in Global system
    slipDirectionG2(iSS,:) = ss(2,:,iSS) * m2;     % for grain 1, [slip direction] expressed in Global system
    slipDirectionG2W(iSS,:) = slipDirectionG2(iSS,:)*burgersWeight(iSS);
    schmidFactorG2(iSS,1) = slipPlaneG2(iSS,:) * stressTensor * slipDirectionG2(iSS,:)';  % Schimid factor for TWIN system j
end


mPrimeMatrix = zeros(nSS,nSS);
resBurgersMatrix = zeros(nSS,nSS);
mPrimeMatrixAbs = zeros(nSS,nSS);
resBurgersMatrixAbs = zeros(nSS,nSS);

sfG1Matrix = repmat(schmidFactorG1,1,nSS);
sfG2Matrix = repmat(schmidFactorG2',nSS,1);
avgSFMatrix = (sfG1Matrix + sfG2Matrix)/2;
for iSSG1 = 1:nSS
    for iSSG2 = 1:nSS
        mPrimeMatrix(iSSG1,iSSG2) = abs(dot(slipPlaneG1(iSSG1,:),slipPlaneG2(iSSG2,:))) * dot(slipDirectionG1(iSSG1,:),slipDirectionG2(iSSG2,:));
        resBurgersMatrix(iSSG1,iSSG2) = norm( slipDirectionG1W(iSSG1,:)-slipDirectionG2W(iSSG2,:) );
        resBurgersMatrixAbs(iSSG1,iSSG2) = min(norm( slipDirectionG1W(iSSG1,:)-slipDirectionG2W(iSSG2,:) ),  norm( slipDirectionG1W(iSSG1,:)+slipDirectionG2W(iSSG2,:) ) );
    end
end
mPrimeMatrixAbs(1:24,1:24) = abs(mPrimeMatrix(1:24,1:24));

