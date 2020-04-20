% SIMULTANEOUS PARAMETER EST. & CONTROL
clc
clear
fs = 2000; % sampling rate
a = 0;
b = 20;
t_array = a:(1/fs):b;
reference_trajectory;
% time-varying control signal
f = 0.5;
U = sin(2*pi*f*t_array);
% plot(t_array,U)
l0=7; l1=20; l2=24;% design parameters

% INITIALIZATION
x0 = [1.16; 0; 0.8976];
w0 = ones(1, 4);
p = [1 1 1 1 1 1]'; true_p = [1 1 1 1 1 1]';
p0 = [34.29; 0;0.1490;0.3341;28.220; 0.2405];
u = 0; flag = 0; track = 0;

% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
half_win_size = 25;
y_dataR = []; y_dot_dataR = []; y_ddot_dataR = []; y_tdot_dataR = [];
y_dataE = []; y_dot_dataE = []; y_ddot_dataE = []; y_tdot_dataE = [];

v=0;
window_num = 1;

y_data = []; y_dot_data = []; y_ddot_data = []; y_tdot_data = [];
x_data = []; U_data = [];
diff = norm(p-p0);
phase = 1;

tic
for i = 26:25:length(t_array)-25
    
    start_pt = i - half_win_size; end_pt = i + half_win_size;
    u = U(start_pt:end_pt);
    [t_arr, x_array] = ode45(@(t,x)sys_sync_gen(t,x,u,t_array(start_pt:end_pt)),t_array(start_pt:end_pt),x0);
    x_array = x_array';
    y = [1 0 0]*x_array;
    y_dT = x_array(2,:);
    y_ddT = gradient(y_dT,(1/fs));
    y_tdT = gradient(y_ddT,(1/fs));
    max_count = 100000;
    y_n = y;
    
    % 4th ORDER LINEAR PARAMETER EST.
    if(i ~= 26)
        w_prev = w_new;
    else
        w_prev = w0;
    end

    ff = @(w)cost_func_sync_gen(w, w_prev, t_array(start_pt: end_pt), y_n, u);
    options = optimoptions('fminunc', 'MaxIterations', 20, 'UseParallel' ,true);
    X = fminunc(ff, w0, options);

    w_new = X;

    a0 = w_new(1);
    a1 = w_new(2);
    a2 = w_new(3);
    a3 = w_new(4);
    disp([a0, a1, a2, a3]);
    estimated_params = [a0,a1,a2,a3];

    % estimation over local window ORDER 4
    disp('Estimating Trajectories with minimized weights');
    [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win, index] = output_estimator(t_array(start_pt:end_pt), y_n, u, t_array(start_pt), t_array(end_pt), estimated_params, max_count);

    y_dataE = [y_dataE y_est_win(index)];
    y_dot_dataE = [y_dot_dataE y_dot_est_win(index)];
    y_ddot_dataE = [y_ddot_dataE y_ddot_est_win(index)];
    y_tdot_dataE = [y_tdot_dataE y_tdot_est_win(index)];
     % capture real data
    y_dataR = [y_dataR y(index)];
    y_dot_dataR = [y_dot_dataR y_dT(index)];
    y_ddot_dataR = [y_ddot_dataR y_ddT(index)];
    y_tdot_dataR = [y_tdot_dataR y_tdT(index)];
    x_data = [x_data t_array(i)];
    U_data = [U_data u(index)];
    x0 = x_array(:,index);
    
    % optimization
    if(diff>0.5)
        [true_p, diff] = optimize_params(p,y_dataE,y_dot_dataE,y_ddot_dataE,y_tdot_dataE,U_data,phase,diff);       
    else
        break;
    end
    p=true_p;
    window_num = window_num + 1
    
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