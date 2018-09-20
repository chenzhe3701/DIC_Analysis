% [abs_schmid_factor, Sorted_SF, burgersXY] = trace_analysis_TiMgAl(varargin)
%
% input can be trace_analysis_Ti(angle sets)
%
% output: 
% abs_schmid_factor = [#, SF, angle_XtoY, trace_x_end, trace_y_end]
% burgersXY = [burgers_X, burgers_Y, ratio].
%
% Zhe Chen, 2015-11-30 revised.
% Zhe Chen, 2016-5-1 revised.  Use parse_angle_inputs().
% Zhe Chen, 2016-6-22 revise to make it for twin.
% Zhe Chen, 2016-6-30, revise to use 1x3 vector for phi_error
%
% chenzhe, 2017-06-07, make it from Ti to TiMgAl.  change input and output

function [abs_schmid_factor, sorted_schmid_factor, burgersXY] = trace_analysis_TiMgAl(euler, phi_sys, phi_error, stressTensor, sampleMaterial, twinTF)

euler = euler + ones(size(euler))*100*eps;
% if TF
% disp(['Euler = [',num2str(euler(1)),',',num2str(euler(2)),',',num2str(euler(3)),'],  phi_sys = ',num2str(phi_sys),',  phi_error = ',num2str(phi_error)]);
% disp([]);

%  define slip systems, format: [slip PLANE, slip DIRECTION]
[ssa, c_a, nss, nts, ~] = define_SS(sampleMaterial,twinTF);
ss = crystal_to_cart_ss(ssa,c_a);

g = euler_to_transformation(euler,phi_sys,phi_error);   % g is the 'transformation matrix' defined in continuum mechanics: x=gX, X is Global corrdinate

for j = 1:1:nss
    N = (g.'*ss(1,:,j).').';  % slip plane (i.e., plane normal) expressed in Global system
    M = (g.'*ss(2,:,j).').';  % slip direction expressed in Global system
    burgers_Z(j) = M*[0 0 1].';    % Z-component of burger's vector of slip system j (Percentage of this vector in Z direction)
    schmid_factor(j,1) = N*stressTensor*M.';  % Schimid factor for slip system j
    abs_schmid_factor(j,1) = j;
    abs_schmid_factor(j,2) = abs(schmid_factor(j,1));  % used to evaluate rank of SF
    sliptrace(j,:) = cross(N, [0 0 1]);
    sliptrace(j,:) = sliptrace(j,:)/norm(sliptrace(j,:));
    trace_X(j,:) = [0, sliptrace(j,1)];  % (trace_X(ii), trace_Y(ii)) are coordinate of end of lines starting from origin representing slip traces
    trace_Y(j,:) = [0, sliptrace(j,2)];
    abs_schmid_factor(j,3) = atand(trace_Y(j,2)/trace_X(j,2));  % angle from X to Y
    abs_schmid_factor(j,4) = trace_X(j,2);                      % trace_end_x
    abs_schmid_factor(j,5) = trace_Y(j,2);                      % trace_end_y
    burgersXY(j,1) = M * [1 0 0]';      % burgers-X component
    burgersXY(j,2) = M * [0 1 0]';      % burgers-Y component
    burgersXY(j,3) = burgersXY(j,1)/burgersXY(j,2);
end
for j = nss+1:nss+nts
    N = (g.'*ss(1,:,j).').';  % slip plane (i.e., plane normal) expressed in Global system
    M = (g.'*ss(2,:,j).').';  % slip direction expressed in Global system
    burgers_Z(j) = M*[0 0 1].';    % Z-component of burger's vector of slip system j (Percentage of this vector in Z direction)
    schmid_factor(j,1) = N*stressTensor*M.';  % Schimid factor for slip system j
    abs_schmid_factor(j,1) = j;
    abs_schmid_factor(j,2) = schmid_factor(j,1);  % used to evaluate rank of SF
    sliptrace(j,:) = cross(N, [0 0 1]);
    sliptrace(j,:) = sliptrace(j,:)/norm(sliptrace(j,:));
    trace_X(j,:) = [0, sliptrace(j,1)];  % (trace_X(ii), trace_Y(ii)) are coordinate of end of lines starting from origin representing slip traces
    trace_Y(j,:) = [0, sliptrace(j,2)];
    abs_schmid_factor(j,3) = atand(trace_Y(j,2)/trace_X(j,2));  % angle from X to Y
    abs_schmid_factor(j,4) = trace_X(j,2);                      % trace_end_x
    abs_schmid_factor(j,5) = trace_Y(j,2);                      % trace_end_y
    burgersXY(j,1) = M * [1 0 0]';      % burgers-X component
    burgersXY(j,2) = M * [0 1 0]';      % burgers-Y component
    burgersXY(j,3) = burgersXY(j,1)/burgersXY(j,2);
end

sorted_schmid_factor = sortrows (abs_schmid_factor, -2);   % sort Schmid Factor: #, SF, angle
% disp([abs_schmid_factor, burgersRU]);

end


