% script to demonstrate the use of my transfer function
%
% chenzhe, 2018-10-03

sf = -0.5:0.01:0.5;
xm = 0.1;
xp = 0.35;
p = 0.9;

% xm correspond to 50% probability
% xp correspond to p probability
yy = transfer_to_logsig(sf, xm, xp, p); 

h1 = figure;
plot(sf, yy);
h2 = figure;
plot(sf, logsig(yy));