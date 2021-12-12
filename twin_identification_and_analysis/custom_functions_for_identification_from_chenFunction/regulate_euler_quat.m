% [q0,q1,q2,q3,phi1New,phiNew,phi2New] = regulate_quaternion(phi1,phi,phi2)
%
% regulate euler-angles and quaternion for HCP
% angles are in degree.

% Zhe Chen, 2015-08-28 revised.

function [q0,q1,q2,q3,phi1New_d,phiNew_d,phi2New_d] = regulate_euler_quat(phi1_d,phi_d,phi2_d)
% phi1=phi1{1};phi=phi{1};phi2=phi2{1};     % for debug
% HCP symmetry operation.  Derived from code by varying EulerAngles.
% To understand, combine with the cubic_symmetry code, and convert to
% angle-axis pairs for better understanding, and better visualization.
disp('This is for HCP');

S = [1, 0, 0, 0;
    sqrt(3)/2, 0, 0, 1/2;
    1/2, 0, 0, sqrt(3)/2;
    0, 0, 0, -1;
    1/2, 0, 0, -sqrt(3)/2;
    sqrt(3)/2, 0, 0, -1/2;
    0, 1, 0, 0;
    0, sqrt(3)/2, -1/2, 0;
    0, 1/2, -sqrt(3)/2, 0;
    0, 0, 1, 0;
    0, 1/2, sqrt(3)/2, 0;
    0, sqrt(3)/2, 1/2, 0];

S = repmat(S,length(phi1_d(:)),1);
quatMat = angle2quat(phi1_d(:)/180*pi,phi_d(:)/180*pi,phi2_d(:)/180*pi,'zxz');
quatMat = kron(quatMat,ones(12,1));

Q = quatmultiply(quatMat,S);

t = 2*(Q(:,4)>0)-1;         % find those with last quaternion element < 0
Q = bsxfun(@times,Q,t);     % make rotation axis pointing upward
index = (Q(:,2)>=0) & (Q(:,3)>=0) & (Q(:,2)>=Q(:,3)*sqrt(3));    % find the one within 0-30 on IPF

index=reshape(index,12,[]);     % if multiple quaternion satisfies criterion (this is rare, but happend once), chose 1st one
index2=sum(index);
indexColumn = find(index2~=1);
for ii = 1:length(indexColumn)
    indC = indexColumn(ii);
    indR = find(index(:,indC)==1);
    index(:,indC) = zeros(12,1);
    index(indR(1),indC) = 1;
end
index = reshape(index,[],1);

q_unique = Q(index,:);                  % find the unique quaternion
q0 = reshape(q_unique(:,1),size(phi1_d));
q1 = reshape(q_unique(:,2),size(phi1_d));
q2 = reshape(q_unique(:,3),size(phi1_d));
q3 = reshape(q_unique(:,4),size(phi1_d));


[phi1New_d,phiNew_d,phi2New_d] = quat2angle([q0(:),q1(:),q2(:),q3(:)],'zxz');     % convert the quaternion into euler angle again
phi1New_d = reshape(phi1New_d, size(q0))/pi*180;
phiNew_d = reshape(phiNew_d, size(q0))/pi*180;
phi2New_d = reshape(phi2New_d, size(q0))/pi*180;

display('regulated euler angle and quaternion');
display(datestr(now));


% % The following checks that the regulated euler angles are the same orientation
% k=round(rand()*700000);
% euler1=[phi1{1}(k),phi{1}(k),phi2{1}(k)]
% euler2=[phi1New(k),phiNew(k),phi2New(k)]
% calculate_misorientation_hcp(euler1,euler2)

% % original code
% for ii=1:length(q0(:))
%     q = [q0(ii),q1(ii),q2(ii),q3(ii)];
%     Q = quatmultiply(q,S);
%     axisAngle = quat2axang(Q);
%     angleAxis = circshift(axisAngle,[0,1]);
%     V = angleAxis(:,2:4);
%     t = 2*(angleAxis(:,4)>0)-1;
%     V = bsxfun(@times,V,t);
%     Q = bsxfun(@times,Q,t);     % make rotation axis pointing upward
%     index = find((V(:,1)>0) & (V(:,2)>0) & (V(:,1)>V(:,2)*sqrt(3)));
%     q_unique = Q(index,:);
%     q0(ii) = q_unique(1);
%     q1(ii) = q_unique(2);
%     q2(ii) = q_unique(3);
%     q3(ii) = q_unique(4);
% end

% % plot to visualize
% figure;hold on;
% for ii=1:12
%     plot3([0,Q(ii,2)],[0,Q(ii,3)],[0,Q(ii,4)])
% end
