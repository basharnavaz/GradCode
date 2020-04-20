% SIMULTANEOUS PARAMETER EST. & CONTROL
% Basharnavaz Khan 
% Earlier version was from Shaunak and Praveen
clc
clear
fs = 1000; % sampling rate
a = 0; b = 10;
t_array = a:(1/fs):b;
% INITIALIZATION
w0 = ones(1, 4);
half_win_size = 25;
u = zeros(1,(half_win_size*2 + 1)); w_new=ones(1,4);
window_num = 1;

y_dataR = []; y_dot_dataR = []; y_ddot_dataR = []; y_tdot_dataR = [];
y_dataE = []; y_dot_dataE = []; y_ddot_dataE = []; y_tdot_dataE = [];
x_data = []; U_data = []; w_matrix = [];

fileID = fopen('y_nl.txt', 'r');
y_trajectory = fscanf(fileID, '%f');
%y_trajectory = awgn(y_trajectory, 20, 'measured');

% Sliding Window 
options = optimoptions('fminunc', 'MaxIterations', 20, 'UseParallel' ,true, 'Display', 'iter-detailed');
for i = 26:25:length(t_array)-25
    
    start_pt = i - half_win_size; end_pt = i + half_win_size;
    t_arr = t_array(start_pt:end_pt);
    y = y_trajectory(start_pt:end_pt);
    y_dT = gradient(y,(1/fs));
    y_ddT = gradient(y_dT,(1/fs));
    y_tdT = gradient(y_ddT,(1/fs));
    max_count = 100000;
    y_n = y;
    %y_n = awgn(y, 40);
    
    % Estimation over local window ORDER 4
    w_prev = w_new;
    ff = @(w)cost_func_sync_gen(w, w_prev, t_arr, y_n, u);
    w_new = fminsearch(ff, w_prev);
    w_matrix = [w_matrix; w_new];
    [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win] = output_point_estimator(t_arr(26), y_n, u, t_arr(1), t_arr(end), w_new, max_count);
    
    % Capture estimated data
    y_dataE = [y_dataE y_est_win];
    y_dot_dataE = [y_dot_dataE y_dot_est_win];
    y_ddot_dataE = [y_ddot_dataE y_ddot_est_win];
    y_tdot_dataE = [y_tdot_dataE y_tdot_est_win];
    
    % Capture real data
    y_dataR = [y_dataR y(26)];
    y_dot_dataR = [y_dot_dataR y_dT(26)];
    y_ddot_dataR = [y_ddot_dataR y_ddT(26)];
    y_tdot_dataR = [y_tdot_dataR y_tdT(26)];
    x_data = [x_data t_array(i)];
    
    window_num = window_num + 1;
    fprintf("Window Number: %d \n", window_num);
end

y_real = interp1(x_data, y_dataR, t_array, 'spline');
y_dot_real = interp1(x_data, y_dot_dataR, t_array, 'spline');
y_ddot_real = interp1(x_data, y_ddot_dataR, t_array, 'spline');
y_tdot_real = interp1(x_data, y_tdot_dataR, t_array, 'spline');

y_est = interp1(x_data, y_dataE, t_array, 'spline');
y_dot_est = interp1(x_data, y_dot_dataE, t_array, 'spline');
y_ddot_est = interp1(x_data, y_ddot_dataE, t_array, 'spline');
y_tdot_est = interp1(x_data, y_tdot_dataE, t_array, 'spline');

plot(t_array,y_est);
