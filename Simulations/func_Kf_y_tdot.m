function [v] = func_Kf_y_tdot(t, t_array, y_array, u_array, t_i, a, a0, a1, a2, a3)
y = interp1(t_array,y_array,t);
u = interp1(t_array,u_array,t);
v = (-24 + 24*a3*(t-a) - 12*a2*(t-a).^2 + 4*a1*(t-a).^3 - a0*(t-a).^4).*y...
    + ((t-a).^4).*u;
end