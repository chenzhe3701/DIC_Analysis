% This is a script to demonstrate the JJ Jonas paper was HKL.
% For me, I should use different equations.

clc;
[ssa, c_a] = define_SS('Mg','twin');
[ss, c_a] = define_SS_cart('Mg','twin');
g = [1 0 0; 0 1 0; 0 0 1];

ssa = ssa(:,:,24)
ss = ss(:,:,24)

n_uvtw = ssa(1,:)
b_uvtw = ssa(2,:)

nu = n_uvtw(1)
nv = n_uvtw(2)
nt = n_uvtw(3)
nw = n_uvtw(4)

bu = b_uvtw(1)
bv = b_uvtw(2)
bt = b_uvtw(3)
bw = b_uvtw(4)


N = ss(1,:) * g;
M = ss(2,:) * g;
MN = M'*N

nU = (2*nu+nv)/sqrt(3)
nV = nv
nW = nw/c_a
bU = (bu+bv/2)*sqrt(3)
bV = 1.5*bv
bW = bw*c_a
qU = nV*bW - nW*bV
qV = nW*bU - nU*bW
qW = nU*bV - nV*bU
disp('-------------')
% but I think should be
nU = nu
nV = (2*nv+nu)/sqrt(3)
nW = nw/c_a
bU = 1.5*bu
bV = (bu/2+bv)*sqrt(3)
bW = bw*c_a
q = cross([nU,nV,nW],[bU,bV,bW])
qU = q(1)
qV = q(2)
qW = q(3)

h = [bU bV bW; qU qV qW; nU nV nW]

h1=h
h1(1,:) = h1(1,:)/norm(h1(1,:));
h1(2,:) = h1(2,:)/norm(h1(2,:));
h1(3,:) = h1(3,:)/norm(h1(3,:));
h1

gamma = 0.129

mine = gamma*MN
his = h1'*[0 0 gamma; 0 0 0; 0 0 0]*h1
his = pinv(h1)*[0 0 gamma; 0 0 0; 0 0 0]*h1
his_wrong = pinv(h)*[0 0 gamma; 0 0 0; 0 0 0]*h
