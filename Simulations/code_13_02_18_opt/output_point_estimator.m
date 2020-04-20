function [y_est, y_dot_est, y_ddot_est, y_tdot_est] = output_point_estimator(t, y_n, u, a, b, estimated_params, max_count)
    % code similar to output_estimator.m 
    % returns the estimated values at the point in the window to be
    % estimated
    
    a0 = estimated_params(1);
    a1 = estimated_params(2);
    a2 = estimated_params(3);
    a3 = estimated_params(4);
    
    t_array = a:(1/1000):b;
    t_i = t;
    n1 = (t_i-a); n2 = (b-t_i);
    com_divisor = n1.^4 + n2.^4;

    % y
    int_K_f = quadgk(@(t,v)myfun_KG_15(t,t_array,y_n,u,t_i, a, a0, a1, a2, a3),a,t_i,'MaxIntervalCount',max_count);
    int_K_b = quadgk(@(t,v)myfun_KG_16(t,t_array,y_n,u,t_i, b, a0, a1, a2, a3),t_i,b,'MaxIntervalCount',max_count);
    int_K = int_K_f + int_K_b;
    y_est = int_K./com_divisor;
    
    % y_dot
    v1 = (-a3*(t_i-a)^4 + 12*(t_i-a)^3) * y_est... 
            + (quadgk(@(t,v)func_Kf_y_dot(t,t_array,y_n,u,t_i,a,a0,a1,a2,a3),a,t_i,'MaxIntervalCount',max_count));
    v2 = (-12*(b-t_i)^3 - a3*(b-t_i)^4) * y_est...
            + (quadgk(@(t,v)func_Kb_y_dot(t,t_array,y_n,u,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count));
    y_dot_est = (v1+v2)./com_divisor;

    % y_ddot
     v3 = (8*(t_i-a).^3 - a3*(t_i-a).^4).*y_dot_est...
         + (-36*(t_i-a).^2 + 8*a3*(t_i-a).^3 - a2*(t_i-a).^4).*y_est...
         + (quadgk(@(t,v)func_Kf_y_ddot(t,t_array,y_n,u,t_i,a,a0,a1,a2,a3),a,t_i,'MaxIntervalCount',max_count));
     v4 = (-8*(b-t_i).^3 - a3*(b-t_i).^4).*y_dot_est...
         + (-36*(b-t_i).^2 - 8*a3*(b-t_i).^3 - a2*(b-t_i).^4).*y_est...
         + (quadgk(@(t,v)func_Kb_y_ddot(t,t_array,y_n,u,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count));
     y_ddot_est = (v3+v4)./com_divisor;

     % y_tdot
     v5 = (4*(t_i-a).^3 - a3*(t_i-a).^4).*y_ddot_est...
         + (-12*(t_i-a).^2 + 4*a3*(t_i-a).^3 - a2*(t_i-a).^4).*y_dot_est...
         + (24*(t_i-a) - 12*a3*(t_i-a).^2 + 4*a2*(t_i-a).^3 - a1*(t_i-a).^4).*y_est...
         + quadgk(@(t,v)func_Kf_y_tdot(t,t_array,y_n,u,t_i,a,a0,a1,a2,a3),a, t_i,'MaxIntervalCount',max_count);
     v6 = (-4*(b-t_i).^3 - a3*(b-t_i).^4).*y_ddot_est...
         + (-12*(b-t_i).^2 - 4*a3*(b-t_i).^3 - a2*(b-t_i).^4).*y_dot_est...
         + (-24*(b-t_i) - 12*a3*(b-t_i).^2 - 4*a2*(b-t_i).^3 - a1*(b-t_i).^4).* y_est...
         + quadgk(@(t,v)func_Kb_y_tdot(t,t_array,y_n,u,t_i,b,a0,a1,a2,a3),t_i, b,'MaxIntervalCount',max_count);
     y_tdot_est = (v5+v6)./com_divisor;

end