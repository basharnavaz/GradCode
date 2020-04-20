function [v] = func_Kb_y_ddot(t, t_array, y_array, u_array, t_i, b, a0, a1, a2, a3)
y = interp1(t_array,y_array,t);
u = interp1(t_array,u_array,t);
v = (96*(b-t) + 36*a3*(b-t).^2 + 8*a2*(b-t).^3 + a1*(b-t).^4 ...
     + (t_i-t).*(24 + 24*a3*(b-t) + 12*a2*(b-t).^2 + 4*a1*(b-t).^3 + a0*(b-t).^4)).*y...
     - (t_i-t).*((b-t).^4).*u;
end