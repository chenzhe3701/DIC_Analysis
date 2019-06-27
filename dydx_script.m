% assume A=1, g(theta) = theta^2

close all;
[r,t] = meshgrid(0.01:0.01:2*pi, 0.01:0.01:2*pi);
% [r,t] = meshgrid(2.3:0.01:2.6, 2.3:0.01:2.6);
y = r.^(9/8).*(t.^2);
x = r.*cos(t);

figure;imagesc(y);
figure;imagesc(x);
A = gradient(y)./gradient(x);
figure;surf(A,'edgecolor','none');

A2 = t.^2 * 9/8.*(r).^(1/8).*(cos(t)).^(-1);
figure;surf(A2,'edgecolor','none');

figure;surf(A-A2,'edgecolor','none');

% A1 = 9/8*r.^(1/8).*(t.^2)./cos(t) - r.^(9/8)*2.*t./r./sqrt(1-cos(t).^2);
% figure;surf(A1,'edgecolor','none');

% figure;surf(A2-A1,'edgecolor','none');


A3 = r.^(1/8)*(9/8.*cos(t).*t.^2 - sin(t).*2.*t);
figure;surf(A3,'edgecolor','none');

