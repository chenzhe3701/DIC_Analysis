% chenzhe, 2018-01-23
% understanding multivariate normal distribution

% generate the position data
[xx,yy]=meshgrid(-1:0.02:1,-1:0.02:1);

% the 'relative' probability, --, the pdf, at these positions
z1=mvnpdf([xx(:),yy(:)]);
z1=reshape(z1,size(xx));
figure;
surf(xx,yy,z1,'edgecolor','none');

% the 'cumulative' probability, --, the cdf, at these positions
z2=mvncdf([xx(:),yy(:)]);
z2=reshape(z2,size(xx));
figure;
surf(xx,yy,z2,'edgecolor','none');

% just add up the pdf values.  Similar shape as the cdf.
z3=cumsum(cumsum(z1,1),2);
figure;
surf(xx,yy,z3,'edgecolor','none');