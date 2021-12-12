
% diffusion problem

c0 = 0;
cs = 1;
x = 50*10^(-9);
D = 1*10^(-18);     % m2/sec
t = 60*60*1;
cxt = c0 + (cs-c0)*(1-erf(x/2/sqrt(D*t)))