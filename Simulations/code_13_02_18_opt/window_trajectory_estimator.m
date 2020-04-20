function [y_est] = window_trajectory_estimator(t_array, y_n, u, a, b, a0, a1, a2, a3, max_count)
    % code similar to nandu_matlab/CH_ReproducingKernels_20 
    dim = length(t_array);
    int_K_f = zeros(dim, 1);
    int_K_b = zeros(dim, 1);
    int_K = zeros(dim, 1);
    y_est = zeros(dim, 1);
    
    for i = 1:length(t_array)

        t_i = t_array(i);
        n1 = (t_i-a); n2 = (b-t_i);
        com_divisor = n1.^4 + n2.^4;
        
        % y
        int_K_f(i) = quadgk(@(t,v)myfun_KG_15(t,t_array,y_n,u,t_i, a, a0, a1, a2, a3),a,t_i,'MaxIntervalCount',max_count);

        int_K_b(i) = quadgk(@(t,v)myfun_KG_16(t,t_array,y_n,u,t_i, b, a0, a1, a2, a3),t_i,b,'MaxIntervalCount',max_count);

        int_K(i) = int_K_f(i) + int_K_b(i);

        y_est(i) = int_K(i)./com_divisor;
    end
end