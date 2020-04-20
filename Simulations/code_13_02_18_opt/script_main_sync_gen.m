clc 
clear
fs = 2000; % sampling rate
a = 0;
b = 1;
t_array = a:(1/fs):b;

% SYNCHRONOUS GENERATOR (NON-LINEAR) MODEL
x0 = [1.16; 0; 0.8976];    % x0=[1; 0; 2] for synchronous generator | x0=[1; 1; 0] for sample system
[t_array, x_array] = ode45(@(t,x)sys_sync_gen(t,x,1),t_array,x0);
x_array = x_array';
y = [1 0 0]*x_array;

% SAMPLE 4TH ORDER MODEL
% x0 = [(0.0523599); (0.0872655); 0; 0]; % for 4th order LTI system
% [t_array, x_array] = ode45(@sys_new_4th_order_LTI_y,t_array,x0);
% x_array = x_array';
% y = [1 0 0 0]*x_array;

y_true = y';
max_count = 100000;
% y_n = awgn(y, 40);
y_n = y;

figure
plot(t_array, y_true);
hold on
y_dot_true = x_array(2,:)';
plot(t_array, y_dot_true);
y_ddot_true = gradient(y_dot_true,(1/fs));
plot(t_array, y_ddot_true);
y_tdot_true = gradient(y_ddot_true,(1/fs));
plot(t_array, y_tdot_true);
% reference_trajectory
% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
half_win_size = 25;
y_est = zeros(1,length(t_array));
y_dot_est = zeros(1,length(t_array));
y_ddot_est = zeros(1,length(t_array));
y_tdot_est = zeros(1,length(t_array));

window_num = 1;
data_pt_array = 26:25:length(t_array)-25;
% y_data = zeros(length(data_pt_array), 1);
% y_dot_data = zeros(length(data_pt_array), 1);
% y_ddot_data = zeros(length(data_pt_array), 1);
% y_tdot_data = zeros(length(data_pt_array), 1);
y_data = [];
y_dot_data = [];
y_ddot_data = [];
y_tdot_data = [];
x_data = []
y_dT = [];
y_dotT = [];
y_ddotT = [];
y_tdotT = [];
tic
w0 = ones(1, 4);
k = 1;
matched_ests = zeros(6,5);
% i = 1;
u = 1;
for i = 26:length(t_array)-25
    
    start_pt = i - half_win_size; end_pt = i + half_win_size;
    
    % parameter estimation over window using ridge regression
    if(i ~= 26)
        w_prev = w_new;
    else
        w_prev = w0;
    end
    
    ff = @(w)cost_func_sync_gen(w, w_prev, t_array(start_pt: end_pt), y_true(start_pt: end_pt)');
%     ff = @(w)cost_func_sync_gen(w, w_prev, t_array(start_pt: end_pt), y_n);
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
    [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win, index] = output_estimator(t_array(start_pt:end_pt), y_n(start_pt:end_pt), t_array(start_pt), t_array(end_pt), estimated_params, max_count);
    %% breaking code
%     if (y_tdot_true(200)-y_tdot_est_win(index) < 10e-5)
%         matched_ests(k,1) = i;
%         matched_ests(k,2) = y_est_win(index);
%         matched_ests(k,3) = y_dot_est_win(index);
%         matched_ests(k,4) = y_ddot_est_win(index);
%         matched_ests(k,5) = y_tdot_est_win(index);
%         k = k+1;
%     end
%     
%     if(k == 9)
%         break;
%     end
    %%
    
    y_est(i) = y_est_win(index);
    y_dot_est(i) = y_dot_est_win(index);
    y_ddot_est(i) = y_ddot_est_win(index);
    y_tdot_est(i) = y_tdot_est_win(index);
    
    y_data = [y_data y_est(i)];
    y_dot_data = [y_dot_data y_dot_est(i)];
    y_ddot_data = [y_ddot_data y_ddot_est(i)];
    y_tdot_data = [y_tdot_data y_tdot_est(i)];
    x_data = [x_data t_array(i)];
%     y_dT = [y_dT y_true];
%     y_dotT = [y_dotT y_dot_true];
%     y_ddotT = [y_ddotT y_ddot_true];
%     y_tdotT = [y_tdotT y_tdot_true];
    window_num = window_num+1
   if(window_num==1200)
       break;
   end
%     disp(window_num);
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
% Envelope
[up,lo] = envelope(y_est,10,'peak');
plot(t_array,up);
plot(t_array,lo);
legend('x1','y-est','up','down');