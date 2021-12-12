
%% simulate Guan paper
[ss,c_a] = define_SS_cart('Mg','Twin');

normals_input = [0 0 1];
for ii = 19:24
    normals_input = [normals_input; ss(1,:,ii)];
end


plot_on_PF(normals_input, [-50 -20 0], 'ND')

%% simulate contraction twin in Jonas Acta 59(2011)2046 paper
close all;

euler_p = [100 75 40];
eulers = euler_p;

for ii = 7:12
    euler_t = euler_by_twin(euler_p, ii, 'Mg');
    eulers = [eulers; euler_t];
end

% This plot is similar to that in the paper
plot_on_PF([0 0 1], eulers, 'ND');
title('reconstruct a figure similar to that in paper','fontweight','normal');

% So, contractin TS #5 is the one chosen in the paper, becuase the
% following plot gives the same position of the basal pole on the PF, as
% shown in the paper
plot_on_PF([0 0 1], eulers(5+1,:), 'ND');  
title('contraction twin system #5 is the active one in paper','fontweight','normal');


% So, we can calculate the theoretical trace direction, and display it
[abs_schmid_factor, sorted_schmid_factor, burgersXY] = ...
    trace_analysis_TiMgAl(eulers(5+1,:),[0 0 0],[0 0 0],[1 0 0; 0 0 0; 0 0 0],'Mg','ctwin'); 

trace_dir = abs_schmid_factor(24+5,3)   % trace direction from x-dir to y-dir, which matches the figure in paper
disp(['the theoretical trace direction of contraction twin system #5 is: ',num2str(trace_dir)]);

% Similarly, if we plot the pole of this contraction twin plan, and use
% eulre angle of the parent grain, it should be similar to the Guan paper
% i.e., connect the plot center and the pole, the line should be
% perpendicular to the twin trace direction.
ss = define_SS_cart('Mg','ctwin');
plot_on_PF(ss(1,:,24+5), eulers(1,:), 'ND');

