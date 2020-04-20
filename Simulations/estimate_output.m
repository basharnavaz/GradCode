function [y_est, y_dot_est, y_ddot_est] = estimate_output(t_array, y_n, a, b, estimated_params, max_count)
    dim = length(t_array);
    int_K_f = zeros(dim, 1);
    int_K_b = zeros(dim, 1);
    int_K = zeros(dim, 1);
    y_est = zeros(dim, 1);
    % Code obtained from /dpg_matlab/CH_ReproducingKernels_24 
    v1 = zeros(dim, 1);
    v2 = zeros(dim, 1);
    y_dot_est = zeros(dim, 1);
    % Code obtained from /kumar_matlab/CH_ReproducingKernels_19
    v3 = zeros(dim, 1);
    v4 = zeros(dim, 1);
    y_ddot_est = zeros(dim, 1);
    
    a0 = estimated_params(1);
    a1 = estimated_params(2);
    a2 = estimated_params(3);
    
    for i = 1:length(t_array)

        t_i = t_array(i);

        int_K_f(i) = quadgk(@(t,v)myfun_KG_25(t,t_array,y_n,t_i, a,a0,a1,a2),a,t_i,'MaxIntervalCount',max_count);

        int_K_b(i) = quadgk(@(t,v)myfun_KG_26(t,t_array,y_n,t_i, b,a0,a1,a2),t_i,b,'MaxIntervalCount',max_count);

        int_K(i) = int_K_f(i) + int_K_b(i);

        y_est(i) = int_K(i)./((t_i-a).^3 + (b-t_i).^3);
        
        % y_dot
        
        v1(i) = 6*(t_i-a).^2.*y_est(i) - a2.*(t_i-a).^3.*y_est(i)...
            + (quadgk(@(t,v)myfun_KG_27(t,t_array,y_n,t_i,a,a0,a1,a2),a,t_i,'MaxIntervalCount',max_count));
        
        v2(i) = -6*(b-t_i).^2.*y_est(i) -a2.*(b-t_i).^3.*y_est(i)...
            + (quadgk(@(t,v)myfun_KG_28(t,t_array,y_n,t_i,b,a0,a1,a2),t_i, b,'MaxIntervalCount',max_count));
        
        y_dot_est(i) = (v1(i)+v2(i))./((t_i-a).^3 + (b-t_i).^3);
        
        % y_ddot
        
         v3(i) = 3*((t_i-a).^2).*y_dot_est(i) - a2.*(t_i-a).^3.*y_dot_est(i)...
             - 6*(t_i-a).*y_est(i) + 3*a2.*(t_i-a).^2.*y_est(i)...
             - a1*((t_i-a).^3).*y_est(i)...
             + (quadgk(@(t,v)myfun_KG_29(t,t_array,y_n,a,a0,a1,a2),a,t_i,'MaxIntervalCount',max_count));
    
         v4(i) = -3*((b-t_i).^2).*y_dot_est(i) - a2.*(b-t_i).^3.*y_dot_est(i)...
             - 6*(b-t_i).*y_est(i) - 3*a2.*(b-t_i).^2.*y_est(i)...
             - a1*((b-t_i).^3).*y_est(i)...
             + (quadgk(@(t,v)myfun_KG_30(t,t_array,y_n,b,a0,a1,a2),t_i, b,'MaxIntervalCount',max_count));
         
         y_ddot_est(i) = (v3(i)+v4(i))./((t_i-a).^3 + (b-t_i).^3);

    end
end