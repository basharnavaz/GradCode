% REFERENCE TRAJECTORY

x_t0 = 1.16; % rad
x_T = 3; x_T = 0.2; x_T = 0.7;% rad

t0 = 0;
T = 20;
fs = 2000;
t_array_u = (t0:1/fs:T)';

r1=252; r2=1050; r3=1800;
r4=1575; r5=700; r6=126;
unit_vector = ones(length(t_array_u), 1);
nabla = (t_array_u - t0)./(T-t0);
phi = nabla.^5 .* (r1.*unit_vector - r2.*nabla + r3.*nabla.^2 -r4.*nabla.^3 + r5.*nabla.^4 -r6.*nabla.^5);
ref_traj = x_t0.*unit_vector + (x_T - x_t0).*phi;

% plot(t_array_u, ref_traj);
% hold on;
ref_traj_dot = gradient(ref_traj,(1/fs));
ref_traj_ddot = gradient(ref_traj_dot,(1/fs));
ref_traj_tdot = gradient(ref_traj_ddot,(1/fs));

% plot(t_array_u, ref_traj_dot);
% plot(t_array_u, ref_traj_ddot);
% plot(t_array_u, ref_traj_tdot);