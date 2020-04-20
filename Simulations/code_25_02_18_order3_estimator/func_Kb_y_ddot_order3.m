function [v] = func_Kb_y_ddot_order3(t, t_array, z_array, t_i, b, l0,l1,l2)
z = interp1(t_array,z_array,t);
v = (6 + 6*l2*(b-t) + 3*l1*(b-t).^2 + l0*(b-t).^3).*z;
end