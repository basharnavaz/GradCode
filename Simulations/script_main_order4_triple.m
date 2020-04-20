clc 

fs = 2000; % sampling rate
a = 0;
b = 2;
t_array = a:(1/fs):b;

% SYNCHRONOUS GENERATOR (NON_LINEAR) MODEL
x0 = [1; 0; 1];    % x0=[1; 0; 2] for synchronous generator | x0=[1; 1; 0] for sample system
[t_array, x_array] = ode45(@sys_sync_gen,t_array,x0);
x_array = x_array';
y = [1 0 0]*x_array;

% SAMPLE 4TH ORDER MODEL
% x0 = [(0.0523599); (0.0872655); 0; 0]; % for 4th order LTI system
% [t_array, x_array] = ode45(@sys_new_4th_order_LTI_y,t_array,x0);
% x_array = x_array';
% y = [1 0 0 0]*x_array;

y = y';
max_count = 100000;
%y_n = awgn(y, 40);
y_n = y;

figure
plot(t_array, x_array(1,:))
hold on
plot(t_array, x_array(2,:))
y_ddot_true = gradient(x_array(2,:),(1/fs));
plot(t_array, y_ddot_true);
y_tdot_true = gradient(y_ddot_true,(1/fs));

% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
win_size_2 = 200;
% win_shift = 50;
% start_pt = 50;
% end_pt = start_pt + win_size;
window_num = 1;
y_est = zeros(length(t_array), 1);
y_dot_est = zeros(length(t_array), 1);
y_ddot_est = zeros(length(t_array), 1);
y_tdot_est = zeros(length(t_array), 1);
% count = 101;
% while (count+win_size_2)<=length(t_array)
tic
% a0=1154.8332;a1=0;a2=156.96;a3=0;
% estimated_params = [a0,a1,a2,a3];
parfor (count = 201:3801)
    % parameter estimation over window
    start_pt = count - win_size_2; end_pt = count + win_size_2;
    [a0,a1,a2,a3] = get_order4_param_estimates(t_array, t_array(start_pt), t_array(end_pt), y_n, max_count);
    disp([a0, a1, a2, a3]);
    estimated_params = [a0,a1,a2,a3];
    % y and y_dot estimation over local window
    [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win] = estimate_order4_output_triple(t_array(start_pt:end_pt), y_n(start_pt:end_pt), t_array(start_pt), t_array(end_pt), estimated_params, max_count);
    y_est(count) = y_est_win(200);
    y_dot_est(count) = y_dot_est_win(200);
    y_ddot_est(count) = y_ddot_est_win(200);
    y_tdot_est(count) = y_tdot_est_win(200);
    window_num = window_num+1
%     disp(window_num);
end
toc
z = plot(t_array, y_est);
zz = plot(t_array, y_dot_est);
zzz = plot(t_array, y_ddot_est);
legend('y-true','y-dot-true','y-ddot-true','y-est','y-dot-est', 'y-ddot-est')
% Envelope
[up,lo] = envelope(y_est,10,'peak');
plot(t_array,up);
plot(t_array,lo);
legend('x1','y-est','up','down');