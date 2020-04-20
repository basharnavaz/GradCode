function [norm_cost] = cost_func_sync_gen(w, w0, t_array, y, u)

    cost = 0;

    a0 = w(1);
    a1 = w(2);
    a2 = w(3);
    a3 = w(4);

    a = t_array(1);
    b = t_array(end);

    lambda = 0.5;
    gamma = 0.25; 

    max_count = 10000;
    [y_est] = window_trajectory_estimator(t_array, y, u, a, b, a0, a1, a2, a3, max_count);
% This part is commented out by Bashar
% if an error occurs uncomment it. 
    if size(y_est)~=size(y)
        y_new=y';
        cost = cost + gamma*norm((w - w0))*norm((w - w0)) + norm(y_new - y_est)*norm(y_new - y_est) + lambda*dot(w, w);
    else 
        cost = cost + gamma*norm((w - w0))*norm((w - w0)) + norm(y - y_est)*norm(y - y_est) + lambda*dot(w, w);
    end
    norm_cost = norm(cost);
    
end

