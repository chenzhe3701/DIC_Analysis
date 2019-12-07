function [exz_iTwin_jMode,exz_iTwin_jSS] = calculate_exz(euler_1, euler_2,sampleMaterial)
% chenzhe, 2019-10-16
% input two euler angles.
% output the expression of 6 twin systems' displacement gradient in the
% 24 slip/twin system -> 5 modes' coordinate in the neighbor grain.
% exz(itwin_in_parent, j_mode_in_neighbor)
%
% Thought: this may also represents how a twin disp grad tensor aligns with
% a slip/twin system in neighbor. (1) Larger means can be mostly
% accommodated by that system.
% (2) 0 means cannot be, or does not require that system, no matter it is
% easy or hard.
% (3) For twin, it can be negative, it means does not align well, cannot
% accommodate well,

g_1 = euler_to_transformation(euler_1,[0,0,0],[0,0,0]);
g_2 = euler_to_transformation(euler_2,[0,0,0],[0,0,0]);
g = g_2 * g_1'; % from grain_1 to grain_2

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);
gamma = 0.1289; % twin shear for Mg

for iTwin = 1:6
    nn = ss(1,:,18+iTwin);
    bb = ss(2,:,18+iTwin);
    yy = cross(nn, bb);
    DCM = [bb; yy; nn];
    dgm_T = [0 0 gamma; 0 0 0; 0 0 0];
    dgm_C = DCM' * dgm_T * DCM;
    dgm_L = g_1' * dgm_C * g_1;
    dgm_c = g_2 * dgm_L * g_2';
    
%     dispGrad = gamma * bb' * nn;
%     dgm = g * dispGrad * g';   % disGrad expressed in grain_2 coordinate
             
    clear exz_iTwin_jMode_kGb;
    for iss = 1:24
        nz = ss(1,:,iss);
        bx = ss(2,:,iss);
        sy = cross(nz, bx);
        
        dcm_crystal_to_ss = [bx;sy;nz];
        
        dgm_ss(:,:,iss) = dcm_crystal_to_ss * dgm_c * dcm_crystal_to_ss';
        exz_thisTwin_jSS(iss) = dgm_ss(1,3,iss);
    end
    % for ss 1-18, take absolute value
    exz_thisTwin_jSS(1:18) = abs(exz_thisTwin_jSS(1:18));
    
    exz_iTwin_jMode(iTwin,1) = max(exz_thisTwin_jSS(1:3));
    exz_iTwin_jMode(iTwin,2) = max(exz_thisTwin_jSS(4:6));
    exz_iTwin_jMode(iTwin,3) = max(exz_thisTwin_jSS(7:12));
    exz_iTwin_jMode(iTwin,4) = max(exz_thisTwin_jSS(13:18));
    exz_iTwin_jMode(iTwin,5) = max(exz_thisTwin_jSS(19:24));
    
    exz_iTwin_jSS(iTwin,:) = exz_thisTwin_jSS;
end


