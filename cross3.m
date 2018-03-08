% chenzhe 2018-03-08, calculate cross product of three 4d vectors
%
% no condition check

function plane = cross3(a,b,c)
    v = ones(4,1);
    A = [v,a(:),b(:),c(:)]';
    B = A\(det(A)*eye(4));
    plane = reshape(B(:,1),1,4);
end