
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
%
% 2019-08-01.  I looked at the code again.
% (1) I previously used a 'weightd' Burgers vector.  Maybe it's not a good
% idea, because slip can be multiples times of unit Burgers vector.  Maybe
% the most important thing is, assume both (cartesian) unit length of slip,
% how big the 'residual' part (difference) need to be.
% (2) The direction. Slip is two direction, so maybe we calculate the resB
% using both of the two directions, and choose the minimum difference.
% Similarly, without considering two directions, in mPrime calculation,
% cosine can be negative.  But we might need to only consider the positive.
% ---> new.  However, if both are twinning, [maybe] we should only consider
% the true direction ??? 
% So, maybe just use the result 'mPrimeMatrixAbs' and 'resBurgersMatrixAbs'.
%
% % [important detail] round to reduce error in do ranking, chenzhe, 2019-10-19
%
% 2019-12-06, just a comment:
% if we want to add weight by the SFs, we might need to:
% normalize m' by (m'-min)/(max-min) = (m'+1)/2
% normalize SFs by (SF+0.5)/1
% then multiply (m'+1)/2 * (SF_a+0.5) * (SF_b+0.5)
% Further, CRSS might also need to be considered.
% These can be done outside of the code.
%
% 2019-12-06  correct error:
% for mPrime between slip and twin systems, first need to
% check relative direction (twin_dir,gb_normal), then make
% (slip_dir,gb_normal) the same relative direction.

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
    slipDirectionG1W(iSS,:) = slipDirectionG1(iSS,:)*burgersWeight(iSS);
    schmidFactorG1(iSS,1) = abs(slipPlaneG1(iSS,:) * stressTensor * slipDirectionG1(iSS,:)');  % Schimid factor for slip system j
end
for iSS = nss+1 : nss+ntwin
    slipPlaneG1(iSS,:) = ss(1,:,iSS) * m1;         % for grain 1, [slip plane normal] expressed in Global system
    slipDirectionG1(iSS,:) = ss(2,:,iSS) * m1;     % for grain 1, [slip direction] expressed in Global system
    slipDirectionG1W(iSS,:) = slipDirectionG1(iSS,:)*burgersWeight(iSS);
    schmidFactorG1(iSS,1) = slipPlaneG1(iSS,:) * stressTensor * slipDirectionG1(iSS,:)';  % Schimid factor for TWIN system j
end



for iSS = 1:nss
    slipPlaneG2(iSS,:) = ss(1,:,iSS) * m2;         % for grain 2, [slip plane normal] expressed in Global system
    slipDirectionG2(iSS,:) = ss(2,:,iSS) * m2;     % for grain 2, [slip direction] expressed in Global system
    slipDirectionG2W(iSS,:) = slipDirectionG2(iSS,:)*burgersWeight(iSS);
    schmidFactorG2(iSS,1) = abs(slipPlaneG2(iSS,:) * stressTensor * slipDirectionG2(iSS,:)');  % Schimid factor for slip system j
end
for iSS = nss+1 : nss+ntwin
    slipPlaneG2(iSS,:) = ss(1,:,iSS) * m2;         % for grain 2, [slip plane normal] expressed in Global system
    slipDirectionG2(iSS,:) = ss(2,:,iSS) * m2;     % for grain 2, [slip direction] expressed in Global system
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
        s1 = slipDirectionG1(iSSG1,:);
        s2 = slipDirectionG2(iSSG2,:);
        n1 = slipPlaneG1(iSSG1,:);
        n2 = slipPlaneG2(iSSG2,:);
        
        % The purpose is to find the best combination of alignment of <dir,dir> and <normal,normal> and see how good they are geometrically compatible:  
        % For slip system, the slip_plan_normal and slip_direction can both have sign changed   
        % For twin system, if change slip_direction, the slip_plan_normal also need to be changed 
        % 
        % Method-1: First check <normal,normal>, then check <dir,dir> 
        % (1) Try to align slip_plan_normals, to make the dot product positive, <normal, normal> < 90 degree.  
        % If change the slip_plan_normal_dir, also need to flip the slip/twin_dir.  
        % (2) Try to align slip/twin_dirs to the same side of the grain boundary, by checking sign(dot(slip,gb_normal))  
        % (2.1) If slip-slip, can change freely
        % (2.2) If slip-twin, first check twin, then change basal freely  
        % (2.3) If twin-twin, we cannot do anything more
        %
        % Method-2: First check <dir,dir>, then check <normal, normal>  
        % (1) Try to align slip/twin_dirs to the same side of the grain boundary, by checking sign(dot(slip,gb_normal))  
        %   Check dir_1 -> check dir_2 and compare -> if different, change     
        %   If change the slip_dir, also need to flip the plan_normal_dir  
        % (2) Try to align slip_plan_normal_dirs, to make the dot product positive, <normal, normal> < 90 degree.  
        % (2.1) If slip-slip, can change freely
        % (2.2) If slip-twin, can change slip freely
        % (2.3) If twin-twin, cannot do anything more
        
        % use method-2
        % (step-1)
        sign_dirAlign_1 = sign(dot(s1,gbNormal));
        sign_dirAlign_2 = sign(dot(s2,gbNormal));
        if (sign_dirAlign_1~=sign_dirAlign_2)
            s2 = -s2;
            n2 = -n2;
        end
        % (step-2)
        sign_normalAlign = sign(dot(n1,n2));
        if (sign_normalAlign==-1)
            if ismember(iSSG1, 1:nss) && ismember(iSSG2, 1:nss)
                n2 = -n2;
            elseif ismember(iSSG1, 1:nss) && ismember(iSSG2, nss+1:nSS)
                n1 = -n1;
            elseif ismember(iSSG1, nss+1:nSS) && ismember(iSSG2, 1:nss) 
                n2 = -n2;
            elseif ismember(iSSG1, nss+1:nSS) && ismember(iSSG2, nss+1:nSS) 
                % nothing we can do
            end
        end
        

        mPrimeMatrix(iSSG1,iSSG2) = dot(n1,n2) * dot(s1,s2);
        resBurgersMatrix(iSSG1,iSSG2) = norm( s1-s2 );
        resBurgersMatrixAbs(iSSG1,iSSG2) = min(norm(s1-s2),  norm(s1+s2) );
        
%         mPrimeMatrix(iSSG1,iSSG2) = abs(dot(slipPlaneG1(iSSG1,:),slipPlaneG2(iSSG2,:))) * dot(slipDirectionG1(iSSG1,:),slipDirectionG2(iSSG2,:));
%         resBurgersMatrix(iSSG1,iSSG2) = norm( slipDirectionG1(iSSG1,:)-slipDirectionG2(iSSG2,:) );
%         resBurgersMatrixAbs(iSSG1,iSSG2) = min(norm( slipDirectionG1(iSSG1,:)-slipDirectionG2(iSSG2,:) ),  norm( slipDirectionG1(iSSG1,:)+slipDirectionG2(iSSG2,:) ) );
    end
end
mPrimeMatrixAbs(1:24,1:24) = abs(mPrimeMatrix(1:24,1:24));

% chenzhe, 2019-08-01.
% Maybe, for twin, don't consider two directions, but let's see. 
resBurgersMatrixAbs(nss+1:nss+ntwin, nss+1:nss+ntwin) = resBurgersMatrix(nss+1:nss+ntwin, nss+1:nss+ntwin);

% weighted metric
weight = sfG1Matrix .* sfG2Matrix;
mPrimeMatrixW = (mPrimeMatrix + 1) .* weight;

% [important detail] round to reduce error in do ranking, chenzhe, 2019-10-19
mPrimeMatrix = round(mPrimeMatrix,4);
resBurgersMatrix = round(resBurgersMatrix,4);
mPrimeMatrixAbs = round(mPrimeMatrixAbs,4);
resBurgersMatrixAbs = round(resBurgersMatrixAbs,4);

