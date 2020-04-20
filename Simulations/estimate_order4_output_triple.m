function [y_est, y_dot_est, y_ddot_est, y_tdot_est] = estimate_order4_output_triple(t_array, y_n, a, b, estimated_params, max_count)
    % code similar to nandu_matlab/CH_ReproducingKernels_20 
    dim = length(t_array);
    int_K_f = zeros(dim, 1);
    int_K_b = zeros(dim, 1);
    int_K = zeros(dim, 1);
    y_est = zeros(dim, 1);
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
    
    a0 = estimated_params(1);
    a1 = estimated_params(2);
    a2 = estimated_params(3);
    a3 = estimated_params(4);
    
    for i = 1:length(t_array)

        t_i = t_array(i);
        n1 = (t_i-a); n2 = (b-t_i);
        com_divisor = n1.^4 + n2.^4;
        
        % y
        int_K_f(i) = quadgk(@(t,v)myfun_KG_15(t,t_array,y_n,t_i, a, a0, a1, a2, a3),a,t_i,'MaxIntervalCount',max_count);

        int_K_b(i) = quadgk(@(t,v)myfun_KG_16(t,t_array,y_n,t_i, b, a0, a1, a2, a3),t_i,b,'MaxIntervalCount',max_count);

        int_K(i) = int_K_f(i) + int_K_b(i);

        y_est(i) = int_K(i)./com_divisor;
        
        % y_dot

        v1(i) = (-a3*(t_i-a)^4 + 12*(t_i-a)^3) * y_est(i)... 
                + (quadgk(@(t,v)func_Kf_y_dot(t,t_array,y_n,t_i,a,a0,a1,a2,a3),a,t_i,'MaxIntervalCount',max_count));
        
        v2(i) = (-12*(b-t_i)^3 - a3*(b-t_i)^4) * y_est(i)...
                + (quadgk(@(t,v)func_Kb_y_dot(t,t_array,y_n,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count));
        
        y_dot_est(i) = (v1(i)+v2(i))./com_divisor;
        
        % y_ddot
        
         v3(i) = (8*(t_i-a).^3 - a3*(t_i-a).^4).*y_dot_est(i)...
             + (-36*(t_i-a).^2 + 8*a3*(t_i-a).^3 - a2*(t_i-a).^4).*y_est(i)...
             + (quadgk(@(t,v)func_Kf_y_ddot(t,t_array,y_n,t_i,a,a0,a1,a2,a3),a,t_i,'MaxIntervalCount',max_count));
    
         v4(i) = (-8*(b-t_i).^3 - a3*(b-t_i).^4).*y_dot_est(i)...
             + (-36*(b-t_i).^2 - 8*a3*(b-t_i).^3 - a2*(b-t_i).^4).*y_est(i)...
             + (quadgk(@(t,v)func_Kb_y_ddot(t,t_array,y_n,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count));
         
         y_ddot_est(i) = (v3(i)+v4(i))./com_divisor;
         
         % y_tdot
         
%          v5(i) = (4*(t_i-a).^3 - a3*(t_i-a).^4 - 4*(b-t_i).^3 - a3*(b-t_i).^4).*y_ddot_est(i) ...
%              + (-12*(t_i-a).^2 + 4*a3*(t_i-a).^3 - a2*(t_i-a).^4 ...
%              - 12*(b-t_i).^2 - 4*a3*(b-t_i).^3 - a2*(b-t_i).^4).*y_dot_est(i) ...
%              + (24*(t_i-a) - 12*a3*(t_i-a).^2 + 4*a2*(t_i-a).^3 - a1*(t_i-a).^4 ...
%              - 24*(b-t_i) - 12*a3*(b-t_i).^2 - 4*a2*(b-t_i).^3 - a1*(b-t_i).*4).* y_est(i) ...
%              + quadgk(@(t,v)func_Kf_y_tdot(t,t_array,y_n,t_i,a,a0,a1,a2,a3),a, t_i,'MaxIntervalCount',max_count) ...
%              + quadgk(@(t,v)func_Kb_y_tdot(t,t_array,y_n,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count);
         
         v5(i) = (4*(t_i-a).^3 - a3*(t_i-a).^4).*y_ddot_est(i)...
             + (-12*(t_i-a).^2 + 4*a3*(t_i-a).^3 - a2*(t_i-a).^4).*y_dot_est(i)...
             + (24*(t_i-a) - 12*a3*(t_i-a).^2 + 4*a2*(t_i-a).^3 - a1*(t_i-a).^4).*y_est(i)...
             + quadgk(@(t,v)func_Kf_y_tdot(t,t_array,y_n,t_i,a,a0,a1,a2,a3),a, t_i,'MaxIntervalCount',max_count);
         
         v6(i) = (-4*(b-t_i).^3 - a3*(b-t_i).^4).*y_ddot_est(i)...
             + (-12*(b-t_i).^2 - 4*a3*(b-t_i).^3 - a2*(b-t_i).^4).*y_dot_est(i)...
             + (-24*(b-t_i) - 12*a3*(b-t_i).^2 - 4*a2*(b-t_i).^3 - a1*(b-t_i).^4).* y_est(i)...
             + quadgk(@(t,v)func_Kb_y_tdot(t,t_array,y_n,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count);
         
         y_tdot_est(i) = (v5(i)+v6(i))./com_divisor;

    end
end