% function g = euler_to_transformation(euler_d[1x3], phi_sys_d[1x3], phi_error[1x3])
% euler is 1x3 row vector, for euler angles, in degrees.
%
% Chenzhe 2016-2-5
%
% Update 2016-02-15
% 
% I come up with a good way to understand the OIM settings.  Now this code
% is modified, to accomodate the phi_sys angle to translate between euler
% angle coordinate and the global coordinate where I can draw the crystal.
%
% TSL software use ZXZ convention for euler angle, and X(A1) is pointing
% down, Y(A2) is pointing right.  This is how the EBSD-detector sees the
% crystal.
%
% The global 'image processing system' is X-right, Y-down, Z-intoPaper.
% However, there is a relationship between the 'imaging system' and the
% 'EBSD system'.
%
% For MSU system, it is OIM setting-1. A1 is along Y-down direction.  A2 is
% along X-right direction. The transformation from 'imaging system' to
% 'EBSD system' can be described by phi_sys = [90,180,0].
% 
% For UM system, the setting was set as the OIM-setting-2.  It is said that
% this is the setting used by 'most modern SEM systems'. A1 is along -Y
% (up) direction.  A2 is along -X (left) direction. There are two ways to
% understand this setting. --> But, if you look at the physical setup, it
% is setting-1: to do EBSD, sample-top on SEM view moved closer.
%
% For UCSB system, the setting is indeed setting-2.  To do EBSD, SEM view
% sample-top moved farther away.
%
% Method-1: global is still the 'image processing system', however, when I
% align the EBSD detector to the global system, the X-of-EBSD detector is
% now pointing up, the Y-of-EBSD detector is pointing left.
%
% Method-2: Suppose I am the EBSD detector, my X is down, Y is right.  When
% I look outside to see the global system ('image processing system'), I
% find that the X is pointing to my left, the Y is pointing to my up.
%
% For OIM-setting-2, phi_sys is [-90, 180, 0].
%
% If global is 'common ploting system', X-right, Y-up, Z-outPaper, then
% phi_sys = [-90, 0, 0], which was what I used before.  But I will not do
% this any more.
%
% update again on 2016-6-30 to take phi_error as 1x3 vector
%
% chenzhe note, 2017-06-07. UM said it's setting-2.  But I think it's
% actually setting-1.  Remember for UM data.


function g = euler_to_transformation(euler_d, phi_sys_d, phi_error)

phi1 = euler_d(1,1);
PHI = euler_d(1,2);
phi2 = euler_d(1,3);
sys1 = phi_sys_d(1,1);
SYS = phi_sys_d(1,2);
sys2 = phi_sys_d(1,3);
err1 = phi_error(1,1);
ERR = phi_error(1,2);
err2 = phi_error(1,3);

% g_error = [
%     cosd(phi_error) cosd(phi_error-90) 0;
%     cosd(phi_error+90) cosd(phi_error) 0;
%     0 0 1];

g_err1 = [
    cosd(err1) cosd(err1-90) 0;
    cosd(err1+90) cosd(err1) 0;
    0 0 1];
g_ERR = [
    1 0 0;
    0 cosd(ERR) cosd(ERR-90);
    0 cosd(ERR+90) cosd(ERR)];
g_err2 = [
    cosd(err2) cosd(err2-90) 0;
    cosd(err2+90) cosd(err2) 0;
    0 0 1];

g_sys1 = [
    cosd(sys1) cosd(sys1-90) 0;
    cosd(sys1+90) cosd(sys1) 0;
    0 0 1];
g_SYS = [
    1 0 0;
    0 cosd(SYS) cosd(SYS-90);
    0 cosd(SYS+90) cosd(SYS)];
g_sys2 = [
    cosd(sys2) cosd(sys2-90) 0;
    cosd(sys2+90) cosd(sys2) 0;
    0 0 1];

g_phi1 = [
    cosd(phi1) cosd(phi1-90) 0;
    cosd(phi1+90) cosd(phi1) 0;
    0 0 1];
g_PHI = [
    1 0 0;
    0 cosd(PHI) cosd(PHI-90);
    0 cosd(PHI+90) cosd(PHI)];
g_phi2 = [
    cosd(phi2) cosd(phi2-90) 0;
    cosd(phi2+90) cosd(phi2) 0;
    0 0 1];
g = g_phi2*g_PHI*g_phi1 * g_sys2*g_SYS*g_sys1 * g_err2 * g_ERR * g_err1;  % g is the 'transformation matrix' defined in continuum mechanics: x=gX, X is Global corrdinate
