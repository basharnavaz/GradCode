clc 
clear
fs = 2000; % sampling rate
a = 0;
b = 1;
t_array = a:(1/fs):b;
matched_ests = zeros(6,6);
max_count = 100000;

w0 = ones(1, 4);
x0 = [1.16; 0; 0.8976];
u = 1;
k = 1;
window_num = 1;
tol = 10e-5;
r1 = -3+6.*rand(400,1);
r2 = -1+2.*rand(400,1);
r3 = -1+2.*rand(400,1);
X_0 = [r1 r2 r3];
i = 1;
% r = 1:20;
while (k<9)
    start_pt = 1; end_pt = 401;
    u = randi(15);
    x0 = X_0(i,:);
    [t_arr, x_array] = ode45(@(t,x)sys_sync_gen(t, x, u), t_array(start_pt:end_pt), x0);
    x_array = x_array';
    y = [1 0 0]*x_array;
    y_true = y';
    y_n = y;
    y_dot_true = x_array(2,:)'; 
    y_ddot_true = gradient(y_dot_true,(1/fs));
    y_tdot_true = gradient(y_ddot_true,(1/fs));
    
    % parameter estimation over window using ridge regression
    
%     ff = @(w)cost_func_sync_gen(w, w_prev, t_array(start_pt: end_pt), y_true(start_pt: end_pt));
    ff = @(w)cost_func_sync_gen(w, w0, t_array(start_pt: end_pt), y_n);
    options = optimoptions('fminunc', 'MaxIterations', 20, 'UseParallel' ,true);
    X = fminunc(ff, w0, options);
    
    w_new = X;
    
    a0 = w_new(1);
    a1 = w_new(2);
    a2 = w_new(3);
    a3 = w_new(4);
    disp([a0, a1, a2, a3]);
    estimated_params = [a0,a1,a2,a3];
    
    % y and y_dot estimation over local window
    disp('Estimating Trajectories with minimized weights');
    [y_est, y_dot_est, y_ddot_est, y_tdot_est, index] = output_estimator(t_array(start_pt:end_pt), y_n, t_array(start_pt), t_array(end_pt), estimated_params, max_count);
    %% breaking code
    if (y_n(200)-y_est(index) < tol &&...
            y_dot_true(200)-y_dot_est(index) < tol &&...
            y_ddot_true(200)-y_ddot_est(index) < tol &&...
            y_tdot_true(200)-y_tdot_est(index) < tol)
        matched_ests(k,2) = y_est(index);
        matched_ests(k,3) = y_dot_est(index);
        matched_ests(k,4) = y_ddot_est(index);
        matched_ests(k,5) = y_tdot_est(index);
        matched_ests(k,6) = u;
%         matched_ests(k,2) = y_n(200);
%         matched_ests(k,3) = y_dot_true(200);
%         matched_ests(k,4) = y_ddot_true(200);
%         matched_ests(k,5) = y_tdot_true(200);    
        disp(k);k = k+1
    end
    i = i+1;
    %%
    window_num = window_num+1
   
%     disp(window_num);
end