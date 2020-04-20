% l2=4.02; % design parameters
% l1=0.0802;
% l0=0.0008;
% 
l0=7; % design parameters
l1=20;
l2=24;

x0 = [1.16; 0; 0.8976];   
[t_array_u, x_array_ctrl] = ode45(@(t,y)controlled_sys_sync_gen(t, y, t_array_u,...
                            ref_traj, ref_traj_dot, ref_traj_ddot, ref_traj_tdot, l0,l1,l2,...
                            b1,b2,b3,b4,E,P),...
                            t_array_u, x0);
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-16);    
C = 340;
[t_array_u, x_array_ctrl] = ode45(@(t,y)controlled_sys_sync_gen(t, y, t_array_u,...
                            ref_traj, ref_traj_dot, ref_traj_ddot, ref_traj_tdot, l0,l1,l2,...
                            C,C,C,C,C,C),...
                            t_array_u, x0);
[dx, ut, t_arr] = controlled_sys_sync_gen([],[],[],[],[],[],[],[],[],[]);
% [dx ut] = controlled_sys_sync_gen(t, y, t_array,...
%                             ref_traj, ref_traj_dot, ref_traj_ddot, ref_traj_tdot, l0,l1,l2);
x_array_ctrl = x_array_ctrl';
y_u = [1 0 0]*x_array_ctrl;

% plotting method

plot(t_array_u,ref_traj)
hold on
plot(t_array_u,y_u')
legend('Ref. trajectory','Controlled output')
grid on
