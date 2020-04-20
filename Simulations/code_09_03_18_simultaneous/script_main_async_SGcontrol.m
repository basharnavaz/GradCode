% CONTROL OF SG

fs = 2000; % sampling rate
a = 0;
b = 20;
t_array = a:(1/fs):b;
reference_trajectory;
% time-varying control signal
f = 3;f1=linspace(1,8,40001);
% f1=fliplr(f1);
% U = sin(2*pi*f*t_array);
Ul = 8*sin(2*pi*f1.*t_array);%.*exp(-0.5*t_array);%+sawtooth(2*pi*20*t_array).*0.2;
% Uh=0.1*cos(2*pi*50*t_array);
U=Ul;
% plot(t_array,U)
l0=7; l1=20; l2=24;% design parameters
alpha = 1;
% INITIALIZATION
x0 = [0; 0; 0.8976];
w0 = ones(1, 4);
p0 = [34.29; 0; 0.1490; 0.3341; 28.220; 0.2405];
true_p=[34.2564; 0; 0.1491;0.3406;27.7363;0.2418];
half_win_size = 25;
u = zeros(1,(half_win_size*2 + 1)); phase = 1; track = 0; u_nl=0; w_new=0;
z_n=0;

window_num = 1;

y_dataR = []; y_dot_dataR = []; y_ddot_dataR = []; y_tdot_dataR = [];
y_dataE = []; y_dot_dataE = []; y_ddot_dataE = []; y_tdot_dataE = [];
x_data = []; U_data = [];

tic
for i = 26:25:length(t_array)-25
    
    start_pt = i - half_win_size; end_pt = i + half_win_size;
    [t_arr, x_array] = ode45(@(t,x)sys_sync_gen(t,x,u,t_array(start_pt:end_pt)),t_array(start_pt:end_pt),x0);
    x_array = x_array';
    y = [1 0 0]*x_array;
    y_dT = x_array(2,:);
    y_ddT = gradient(y_dT,(1/fs));
    y_tdT = gradient(y_ddT,(1/fs));
    max_count = 100000;
    y_n = y;
%     y_n = awgn(y, 40);


    z_n = y_n - ref_traj(start_pt:end_pt)';
    [y_est_win, y_dot_est_win, y_ddot_est_win, y_tdot_est_win, index] = output_estimator_LE(t_array(start_pt:end_pt), z_n, t_array(start_pt), t_array(end_pt), l0, l1, l2,...
        ref_traj(start_pt:end_pt), ref_traj_dot(start_pt:end_pt), ref_traj_ddot(start_pt:end_pt), ref_traj_tdot(start_pt:end_pt), max_count);

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
    
    F = y_est_win(index);
    F_dot = y_dot_est_win(index);
    F_ddot = y_ddot_est_win(index);
    F_tdot = y_tdot_est_win(index);
    x0 = x_array(:,index);
    
    % compute control
    r = ref_traj(i);
    r_dot = ref_traj_dot(i);
    r_ddot = ref_traj_ddot(i);
    r_tdot = ref_traj_tdot(i);
    v = r_tdot - l2*(F_ddot - r_ddot)- l1*(F_dot - r_dot) - l0*(F - r);
    b1 = true_p(1); b2 = true_p(2); b3 = true_p(3); b4 = true_p(4); P = true_p(5); E = true_p(6);
    u_nl = E - b4*(P - b2*F_dot - F_ddot)/(b1*sin(F)) + b3*cos(F)...
        - ((-b2*F_ddot - v)*sin(F)-(P-b2*F_dot-F_ddot)*F_dot*cos(F))/(b1*(sin(F))^2);
    
    % include saturation HERE if reqd
    
    u_nl_vec = ones(1,(half_win_size*2) + 1)*u_nl;
    u_i_vec = ones(1,(half_win_size*2) + 1)*U(i);
    
 
    u = u_nl_vec;
    window_num = window_num + 1
    
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
hold on;
plot(t_array,ref_traj);
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
grid on
legend({'Controlled output $\mathbf{y(t)}$', 'Reference Trajectory $\mathbf{F^{*}(t)}$'},'Interpreter','latex');

uu = interp1(x_data, U_data, t_array, 'spline');
plot(t_array, uu);
grid on;
legend({'Control signal $\mathbf{u^{*}(t)}$'},'Interpreter','latex');
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
% % plotting code for y and its derivatives
% plot(t_array,y_dr);plot(t_array,y_dest);
% plot(t_array,y_dr);hold on;plot(t_array,y_dest);
% plot(t_array,y_ddr);hold on;plot(t_array,y_ddest);
% plot(t_array,y_tdr);hold on;plot(t_array,y_tdest);
% % extrapolation
% y_ext = interp1(x_data, y_data, t_array(1: 201), 'spline','extrap');
% y_dot_ext = interp1(x_data, y_dot_data, t_array(1: 201), 'spline','extrap');
% y_ddot_ext = interp1(x_data, y_ddot_data, t_array(1: 201), 'spline','extrap');
% y_tdot_ext = interp1(x_data, y_tdot_data, t_array(1: 201), 'spline','extrap');