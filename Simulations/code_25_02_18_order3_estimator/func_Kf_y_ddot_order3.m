function [v] = func_Kf_y_ddot_order3(t, t_array, z_array, t_i, a, l0,l1,l2)
z = interp1(t_array,z_array,t);
v = (6 - 6*l2*(t-a) + 3*l1*(t-a).^2 - l0*(t-a).^3).*z;
end