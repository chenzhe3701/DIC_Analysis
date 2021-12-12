% Chenzhe 2016-2-5
% function ss = hex_to_cart_ss(ssa)
% size(ssa) = 2, 4, 24 or other, defining ss using four indices
% row_1 = plane_normal, row_2 = slip direction, page = slip system
%
% chenzhe, 2017-06-08
% modified from hex_to_cart_ss, make it work for cubic.

function ss = crystal_to_cart_ss(ssa,c_a)
if size(ssa,2) == 4
    for i = 1:size(ssa,3)            % Change n & m to unit vector
        n = [ssa(1,1,i)  (ssa(1,2,i)*2+ssa(1,1,i))/sqrt(3)  ssa(1,4,i)/c_a]; % Plane normal /c_a, into a Cartesian coordinate
        m = [ssa(2,1,i)*3/2  (ssa(2,1,i)+ssa(2,2,i)*2)*sqrt(3)/2  ssa(2,4,i)*c_a]; % Slip direction *c_a, into a Cartesian coordinate
        ss(1,:,i) = n/norm(n);  % slip plane, normalized
        ss(2,:,i) = m/norm(m);  % slip direction, normalized
    end
elseif size(ssa,2) ==3
    for i = 1:size(ssa,3)
        n = ssa(1,:,i);
        m = ssa(2,:,i);
        ss(1,:,i) = n/norm(n);  % slip plane, normalized
        ss(2,:,i) = m/norm(m);  % slip direction, normalized
    end
end

end