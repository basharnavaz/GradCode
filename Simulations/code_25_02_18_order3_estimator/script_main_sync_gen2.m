
fs = 2000; % sampling rate
a = 0;
b = 20;
t_array = a:(1/fs):b;

% SYNCHRONOUS GENERATOR (NON-LINEAR) MODEL
x0 = [1.16; 0; 0.8976];    % x0=[1; 0; 2] for synchronous generator | x0=[1; 1; 0] for sample system

max_count = 100000;

% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
half_win_size = 25;
y_est = zeros(length(t_array), 1);
y_dot_est = zeros(length(t_array), 1);
y_ddot_est = zeros(length(t_array), 1);
y_tdot_est = zeros(length(t_array), 1);

window_num = 1;

y_data = [];
y_dot_data = [];
y_ddot_data = [];
y_tdot_data = [];
x_data = [];
u_data = [];
tic
% b1 = 34.29; b2=0; b3=0.1490; b4=0.3341; P=28.220; E=0.2405;
b1 = 33.1989; b2=0; b3=0.1491; b4=0.3434; P=28.2149; E=0.2421;
b1 = 34.2135; b2=0.0014; b3=0.1577; b4=0.3367; P=28.2293; E=0.2399;
w0 = ones(1, 4);
k = 1;
matched_ests = zeros(6,5);
u = zeros(1,51); u_array = zeros(length(t_array), 1);
l0=7; % design parameters
l1=20;
l2=24;
for i = 26:25:length(t_array)-25
    
    start_pt = i - 25; end_pt = i + 25;
    [t_arr, x_array] = ode45(@(t,x)sys_sync_gen(t,x,u,t_array(start_pt:end_pt)), t_array(start_pt:end_pt), x0);
    x_array = x_array';
    y = [1 0 0]*x_array;
    y_n = y;
    y_dT = x_array(2,:);
    y_ddT = gradient(y_dT,(1/fs));
    y_tdT = gradient(y_ddT,(1/fs));
    if(i<-1)
%         parameter estimation over window using ridge regression
        if(i ~= 201)
            w_prev = w_new;
        else
            w_prev = w0;
        end
        
        ff = @(w)cost_func_sync_gen(w, w_prev, t_array(start_pt: end_pt), y_n);
        options = optimoptions('fminunc', 'MaxIterations', 20, 'UseParallel', true);
        X = fminunc(ff, w0, options);
        
        w_new = X;
        
        a0 = w_new(1);
        a1 = w_new(2);
        a2 = w_new(3);
        a3 = w_new(4);
        disp([a0, a1, a2, a3]);
        estimated_params = [a0,a1,a2,a3];
        
        % estimation over local window
        disp('Estimating Trajectories with minimized weights');
        [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win, index] = output_estimator(t_array(start_pt:end_pt), y_n, t_array(start_pt), t_array(end_pt), estimated_params, max_count);
        %     y_n is the measurement == F
    else
        z_n = y_n - ref_traj(start_pt:end_pt)';
        [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win, index] = output_estimator_LE(t_array(start_pt:end_pt), z_n, t_array(start_pt), t_array(end_pt), l0, l1, l2,...
            ref_traj(start_pt:end_pt), ref_traj_dot(start_pt:end_pt), ref_traj_ddot(start_pt:end_pt), max_count);
    end
    x0 = x_array(:, index);
    %%
    y_est(i) = y_est_win(index); F = y_est_win(index);
    y_dot_est(i) = y_dot_est_win(index); F_dot = y_dot_est_win(index);
    y_ddot_est(i) = y_ddot_est_win(index); F_ddot = y_ddot_est_win(index);
    y_tdot_est(i) = y_tdot_est_win(index); F_tdot = y_tdot_est_win(index);
    
    y_data = [y_data y_est(i)];
    y_dot_data = [y_dot_data y_dot_est(i)];
    y_ddot_data = [y_ddot_data y_ddot_est(i)];
    y_tdot_data = [y_tdot_data y_tdot_est(i)];
    x_data = [x_data t_array(i)];
    u_data = [u_data u];
    window_num = window_num+1
    
    %  compute control signal
    r = ref_traj(i);
    r_dot = ref_traj_dot(i);
    r_ddot = ref_traj_ddot(i);
    r_tdot = ref_traj_tdot(i);
    
    v = r_tdot - l2*(F_ddot - r_ddot)...
        - l1*(F_dot - r_dot) - l0*(F - r);
    
    u = E - b4*(P - b2*F_dot - F_ddot)/(b1*sin(F)) + b3*cos(F)...
        - ((-b2*F_ddot - v)*sin(F)-(P-b2*F_dot-F_ddot)*F_dot*cos(F))/(b1*(sin(F))^2)...
%         +(r-y_est(i));
    u = ones(1,(half_win_size*2) + 1)*u;
%     pause(0.3)
end
y_est = interp1(x_data, y_data, t_array(201:end), 'spline');
y_dot_est = interp1(x_data, y_dot_data, t_array(201:end), 'spline');
y_ddot_est = interp1(x_data, y_ddot_data, t_array(201:end), 'spline');
y_tdot_est = interp1(x_data, y_tdot_data, t_array(201:end), 'spline');

for i = 1: length(y_est)
    if(isnan(y_est(i)))
        y_est(i) = 0;
        y_dot_est(i) = 0;
        y_ddot_est(i) = 0;
        y_tdot_est(i) = 0;
    end
end

% extrapolation
y_ext = interp1(x_data, y_data, t_array(1: 201), 'spline','extrap');
y_dot_ext = interp1(x_data, y_dot_data, t_array(1: 201), 'spline','extrap');
y_ddot_ext = interp1(x_data, y_ddot_data, t_array(1: 201), 'spline','extrap');
y_tdot_ext = interp1(x_data, y_tdot_data, t_array(1: 201), 'spline','extrap');

y_est = [zeros(200,1); y_est];
y_dot_est = [zeros(200,1); y_dot_est];
y_ddot_est = [zeros(200,1); y_ddot_est];
y_tdot_est = [zeros(200,1); y_tdot_est];

for i = 1: 201
    y_est(i) = y_ext(i);
    y_dot_est(i) = y_dot_ext(i);
    y_ddot_est(i) = y_ddot_ext(i);
    y_tdot_est(i) = y_tdot_ext(i);
end
toc
z = plot(t_array, y_est);
zz = plot(t_array, y_dot_est);
zzz = plot(t_array, y_ddot_data);
legend('y-true','y-dot-true','y-ddot-true','y-est','y-dot-est', 'y-ddot-est')
