function euler = random_euler(method)
% Randomly generate an orientation, and output the euler angle.
% % I think method 1 is better than method0, but double check.
%
% chenzhe, 2019-11-18

if exist('method','var')&&(method==1)
    % This is bad by simply put 3 euler angles.
%     phi1 = -180 + rand(1)*360;   % [-pi, pi]
%     phi = rand(1)*180;  % [0 pi]
%     phi2 = -180 + rand(1)*360;
    
elseif exist('method','var')&&(method==2)
    % random rotation about random axis.
    % This seems bad.   Better to define cAxis rather than rotation axis.  
%     z = -1 + rand(1)*2;
%     a = sqrt(1^2-z^2);
%     rho = rand(1) * 360;
%     x = a*cosd(rho);
%     y = a*sind(rho);
%     
%     rotAxis = [x,y,z];
%     rotAngle = rand(1) * 2*pi;
%     
%     rod(1) = tan(1/2*rotAngle)*x;
%     rod(2) = tan(1/2*rotAngle)*y;
%     rod(3) = tan(1/2*rotAngle)*z;
%     
%     [R1, R2, R3] = rod2angle(rod,'zxz');
%     phi1 = R1/pi*180;
%     phi = R2/pi*180;
%     phi2 = R3/pi*180;
    
else
    % Directly define the 3 axis' direction
    %
    % spherical cap, area = 2pi*r*h = pi(a^2+h^2) = 2pi*r^2*(1-cos(t))
    % r = sphere radius, h = height of cap, a = base radius
    % t = polar angle between rays from the center of the sphere to the apxe of
    % the cap (pole) and the edge of the disk forming the based of the cap
    %
    % let z = r-h, the area of the ring is 2pi*r*(r-h) = 2pi*rz.
    % Therefore, it's random along z parameter.
    % Therefore, z=rand(1), rotation angle = rand(1)*360, can describe a random orientation
    
    z = -1 + rand(1)*2;
    a = sqrt(1^2-z^2);
    rho = rand(1) * 360;
    x = a*cosd(rho);
    y = a*sind(rho);
    
    cAxis = [x,y,z];    % random cAxis of crystal -- get this by [phi1, phi]  
    phi1 = atand(-x/y);
    phi = atand(sqrt(x^2+y^2)/z);
    
    Angle = rand(1) * 360;  % random rotation of crystal about cAxis
    phi2 = Angle;

end

euler = [phi1,phi,phi2];
% xyz = [x,y,z];

% can use this to check.  Note that pole figure is 'equal area' projection   
if 0
    close all;
    clear euler euler_1 euler_2
    for ii=1:100
%         euler_1(ii,:) = random_euler(1);
%         euler_2(ii,:) = random_euler(2);
        euler(ii,:) = random_euler();
    end
%     plot_on_PF([0 0 1], euler_1, 'ND');
%     plot_on_PF([0 0 1], euler_2, 'ND');
    plot_on_PF([0 0 1], euler, 'ND');
%     figure;boxplot(euler_1);
%     figure;boxplot(euler_2);
%     figure;boxplot(euler);
end

% randomly generate 2 orientation, cross to get the 3rd, to get dcm.
% This is bad, because cannot garantee x and y orthogonal.
%     dcm = zeros(3);
%
%     z = -1 + rand(1)*2;
%     a = sqrt(1^2-z^2);
%     rho = rand(1) * 360;
%     x = a*cosd(rho);
%     y = a*sind(rho);
%
%     dcm(1,:) = [x,y,z];
%
%     z = -1 + rand(1)*2;
%     a = sqrt(1^2-z^2);
%     rho = rand(1) * 360;
%     x = a*cosd(rho);
%     y = a*sind(rho);
%
%     dcm(2,:) = [x,y,z];
%     dcm(3,:) = cross(dcm(1,:),dcm(2,:));
%
%     [R1, R2, R3] = dcm2angle(dcm,'zxz');
%     phi1 = R1/pi*180;
%     phi = R2/pi*180;
%     phi2 = R3/pi*180;