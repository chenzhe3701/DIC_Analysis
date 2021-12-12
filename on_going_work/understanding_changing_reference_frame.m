[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
ss = crystal_to_cart_ss(ssa,c_a);
gamma = 0.1289; % twin shear for Mg

euler = [-178.0300   57.1230  142.3200];
g = euler_to_transformation(euler,[0,0,0],[0,0,0]);

n = ss(1,:,19)
b = ss(2,:,19)

N = n * g
B = b * g

f = eye(3) + gamma * b' * n
e = (f'*f-eye(3))/2

F = eye(3) + gamma * B' * N
E = (F'*F-eye(3))/2

g * E * g'  % ok
g' * E * g