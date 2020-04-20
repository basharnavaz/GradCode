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
% U = sin(2*pi*f*t_array);
Ul = sin(2*pi*f*t_array);%.*exp(-0.5*t_array);%+sawtooth(2*pi*20*t_array).*0.2;
% Uh=0.1*cos(2*pi*50*t_array);
U=Ul;
% plot(t_array,U)
l0=7; l1=20; l2=24;% design parameters
alpha = 0;
% INITIALIZATION
x0 = [1.16; 0; 0.8976];
w0 = ones(1, 4);
p = [1 1 1 1 1 1]'; true_p = [1 1 1 1 1 1]';
p0 = [34.29; 0;0.1490;0.3341;28.220; 0.2405];
u = 0; phase = 1; track = 0; u_nl=0;

% PARAMETER AND OUTPUT ESTIMATION OVER DRAGGING WINDOW SNIPPET
half_win_size = 25;
y_est = zeros(1,length(t_array));
y_dot_est = zeros(1,length(t_array));
y_ddot_est = zeros(1,length(t_array));
y_tdot_est = zeros(1,length(t_array));
window_num = 1;
F = 1;
F_dot = 1;
F_ddot = 1;
F_tdot = 1;

y_data = []; y_dot_data = []; y_ddot_data = []; y_tdot_data = [];
x_data = []; U_data = [];
diff = norm(p-p0);

tic
for i = 26:25:length(t_array)-25
    
    start_pt = i - half_win_size; end_pt = i + half_win_size;
    
    % compute control 
    
    r = ref_traj(i);
    r_dot = ref_traj_dot(i);
    r_ddot = ref_traj_ddot(i);
    r_tdot = ref_traj_tdot(i);
    v = r_tdot - l2*(F_ddot - r_ddot)- l1*(F_dot - r_dot) - l0*(F - r);
    b1 = true_p(1); b2 = true_p(2); b3 = true_p(3); b4 = true_p(4); P = true_p(5); E = true_p(6);
    u_nl = E - b4*(P - b2*F_dot - F_ddot)/(b1*sin(F)) + b3*cos(F)...
        - ((-b2*F_ddot - v)*sin(F)-(P-b2*F_dot-F_ddot)*F_dot*cos(F))/(b1*(sin(F))^2);
   
%     if(u_nl>50)
%         u_nl=50;
%     elseif(u_nl<-50)
%         u_nl=-50;
%     end
    u_nl_vec = ones(1,(half_win_size*2) + 1)*u_nl;
    u_i_vec = ones(1,(half_win_size*2) + 1)*U(i);
    switch phase
        case 1
            u = u_i_vec + u_nl_vec ;
%             alpha = alpha-0.005;
        case 2
%             u = exp(-alpha)*u_i_vec + u_nl_vec;
%             u = u_nl_vec;
            u = u_nl_vec;
            alpha = alpha+0.1;
            p=true_p;
        case 3
            u = ones(1,(half_win_size*2) + 1)*u_nl;
            p = true_p;
    end
    
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
    
    y_data = [y_data y_est_win(index)];
    y_dot_data = [y_dot_data y_dot_est_win(index)];
    y_ddot_data = [y_ddot_data y_ddot_est_win(index)];
    y_tdot_data = [y_tdot_data y_tdot_est_win(index)];
    x_data = [x_data t_array(i)];
    U_data = [U_data u(index)];
    F = y_est_win(index);
    F_dot = y_dot_est_win(index);
    F_ddot = y_ddot_est_win(index);
    F_tdot = y_tdot_est_win(index);
    x0 = x_array(:,index);
    
    % optimization
    if(window_num <= 800)
        opt_func = @(p)opt_func_SG_params(p, y_data, y_dot_data, y_ddot_data,...
            y_tdot_data, U_data, phase);
        options = optimoptions('fmincon', 'MaxIterations', 10000,...
            'MaxFunEvals', 10000, 'UseParallel', true, 'Display', 'iter-detailed',...
            'Algorithm', 'sqp');
        A=[]; b=[]; Aeq=[]; beq=[]; nonlcon = [];
        lb = zeros(6,1);
        ub=[100;100;100;100;100;100];
        
        [p,fval] = fmincon(opt_func,p,A,b,Aeq,beq,lb,ub,nonlcon,options);
        
        diff_new = norm(p-p0);
        if(diff_new < diff)
            true_p = p;
            diff = diff_new;
            if(diff_new < 1)
                phase = 2; % estimates quite close to actual
            end
        end
%     else
%         phase = 3;
    end
   
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