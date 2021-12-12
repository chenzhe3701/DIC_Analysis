
%% script to artificially create legends
hold on;
legendString = {'iE=2-3','iE=3-4','iE=4-5'};
N = length(legendString);
colors = lines(N);

% h = zeros(N,1);
for ii = 1:N
    h(ii) = plot(nan,nan,'.','color', colors(ii,:));
end

legend(h,legendString);