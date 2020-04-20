clc 

fs = 1000; % sampling rate
a = 0;
b = 0.5;
x0 = [1; 0; 2];    % x0=[1; 0; 1] for synchronous generator | x0=[1; 1; 0] for sample system
t_array = a:(1/fs):b;
% SYNCHRONOUS GENERATOR (NON_LINEAR) MODEL
[t_array, x_array] = ode45(@sys_sync_gen,t_array,x0);

% SAMPLE 3RD ORDER MODEL
% [t_array, x_array] = ode45(@sys_ThirdOrder,t_array,x0);
x_array = x_array';
y = [1 0 0]*x_array;
y = y';
max_count = 100000;
% y_n = awgn(y, 30);
y_n = y;

figure
plot(t_array, x_array(1,:))
hold on
plot(t_array, x_array(2,:))
y_ddot_true = gradient(x_array(2,:),(1/fs));
plot(t_array, y_ddot_true);
% plot(t_array, x_array(2,:))
% plot(t_array, x_array(3,:))

% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
win_size_2 = 100; % needs to be an odd number
% win_shift = 50;
% start_pt = 50;
% end_pt = start_pt + win_size;
window_num = 1;
y_est = zeros(length(t_array), 1);
y_dot_est = zeros(length(t_array), 1);
y_ddot_est = zeros(length(t_array), 1);
% count = 101;
% while (count+win_size_2)<=length(t_array)
tic
parfor (count = 101:401)
    % parameter estimation over window
    start_pt = count - win_size_2; end_pt = count + win_size_2;
    [a0,a1,a2] = get_param_estimates(t_array, t_array(start_pt), t_array(end_pt), y_n, max_count);
    disp([a0, a1, a2]);
    estimated_params = [a0,a1,a2];
    % FILE WRITE OPERATION
    %     dlmwrite('window_estimates.txt', [a0, a1, a2], '-append', 'delimiter', ' ');
    
    % y and y_dot estimation over local window
    [y_est_win, y_dot_est_win, y_ddot_est_win] = estimate_output(t_array(start_pt:end_pt), y_n(start_pt:end_pt), t_array(start_pt), t_array(end_pt), estimated_params, max_count);
    y_est(count) = y_est_win(100);
    y_dot_est(count) = y_dot_est_win(100);
    y_ddot_est(count) = y_ddot_est_win(100);
    
%     % output estimation over full interval
%     y_est_window = estimate_output(t_array, y_n, a, b, estimated_params, max_count);
%     y_est(start_pt:start_pt + win_shift) = y_est_window(start_pt:start_pt + win_shift);
    
    % loop updates
%     start_pt = start_pt + win_shift;
%     end_pt = end_pt + win_shift;
%     count = count + 1;
%     if(end_pt > length(t_array))
%         break;
%     end
    window_num = window_num+1
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