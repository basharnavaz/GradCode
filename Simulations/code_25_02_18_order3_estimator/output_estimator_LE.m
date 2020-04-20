function [y_est, y_dot_est, y_ddot_est, y_tdot_est, index] = output_estimator_LE(t_array, z_n, a, b, l0,l1,l2,R,R_dot,R_ddot,R_tdot, max_count)

dim = length(t_array);
int_K_f = zeros(dim, 1);
int_K_b = zeros(dim, 1);
int_K = zeros(dim, 1);
y_est = zeros(dim, 1);
summation = zeros(dim, 1);
%
v1 = zeros(dim, 1);
v2 = zeros(dim, 1);
y_dot_est = zeros(dim, 1);
%
v3 = zeros(dim, 1);
v4 = zeros(dim, 1);
y_ddot_est = zeros(dim, 1);
%
v5 = zeros(dim, 1);
v6 = zeros(dim, 1);
y_tdot_est = zeros(dim, 1);

parfor i = 1:length(t_array)
    
    t_i = t_array(i);
    n1 = (t_i-a); n2 = (b-t_i);
    com_divisor = n1.^3 + n2.^3;
    
    % y
    
    int_K_f(i) = quadgk(@(t,v)func_Kf_y_order3(t,t_array,z_n,t_i, a, l0,l1,l2),a,t_i,'MaxIntervalCount',max_count);
    
    int_K_b(i) = quadgk(@(t,v)func_Kb_y_order3(t,t_array,z_n,t_i, b, l0,l1,l2),t_i,b,'MaxIntervalCount',max_count);
    
    int_K(i) = (com_divisor*R(i)) + int_K_f(i) + int_K_b(i);
    
    y_est(i) = int_K(i)./com_divisor;
    
    % COMPUTATION OF KERNEL SQUARED
    forward_K_squared = quadgk(@(t, f)squared_Kf_y_order3(t, t_i, a, l0, l1, l2), a, t_i,'MaxIntervalCount', max_count);
    backward_K_squared = quadgk(@(t, f)squared_Kb_y_order3(t, t_i, b, l0, l1, l2), t_i, b,'MaxIntervalCount', max_count);
    
    summation(i) = (forward_K_squared + backward_K_squared)./com_divisor;
    
    % y_dot
    z_est = y_est(i)-R(i);
    
    v1(i) = (quadgk(@(t,v)func_Kf_y_dot_order3(t,t_array,z_n,t_i,a,l0,l1,l2),a,t_i,'MaxIntervalCount',max_count));
    
    v2(i) = (quadgk(@(t,v)func_Kb_y_dot_order3(t,t_array,z_n,t_i,b,l0,l1,l2),t_i, b,'MaxIntervalCount',max_count));
    
    y_dot_est(i) = (v1(i) + v2(i) +(-l2*com_divisor + 6*((t_i-a).^2 - (b-t_i).^2)).*z_est...
        + com_divisor*R_dot(i))./com_divisor;
    
    % y_ddot
    z_dot_est = y_dot_est(i)-R_dot(i);
    v3(i) = (quadgk(@(t,v)func_Kf_y_ddot_order3(t,t_array,z_n,t_i,a,l0,l1,l2),a,t_i,'MaxIntervalCount',max_count));
    
    v4(i) = (quadgk(@(t,v)func_Kb_y_ddot_order3(t,t_array,z_n,t_i,b,l0,l1,l2),t_i, b,'MaxIntervalCount',max_count));
    
    y_ddot_est(i) = (v3(i) + v4(i) + (3*(t_i-a)^2 - 3*(b-t_i)^2 - l2*com_divisor)*z_dot_est...
        + (-6*(t_i-a) + 3*l2*(t_i-a)^2 - l1*(t_i-a)^3 - 6*(b-t_i) - 3*l2*(b-t_i)^2 - l1*(b-t_i)^3)*z_est...
        + com_divisor*R_ddot(i))./com_divisor;
    
    % y_tdot
    
    y_tdot_est(i) = -(l2*y_ddot_est(i) + l1*y_dot_est(i) + l0*y_est(i))...
                    +(R_tdot(i) + l2*R_ddot(i) + l1*R_dot(i) + l0*R(i));
    
end
[~, index] = min(summation);
end