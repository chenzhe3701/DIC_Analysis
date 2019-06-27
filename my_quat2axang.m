function axang = my_quat2axang(quat)
% the function quat2axang is no longer available.  Therefore, write my own.
% quat = a + bi + cj + dk
% chenzhe, 2019-04-02

g=quat2dcm(quat);
[thetad, v] = derotation(g);

axang(1:3) = v(1:3);
axang(4) = thetad/180*pi;

end
