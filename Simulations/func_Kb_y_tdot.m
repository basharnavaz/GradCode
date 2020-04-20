function [v] = func_Kb_y_tdot(t, t_array, y_array, u_array, t_i, b, a0, a1, a2, a3)
y = interp1(t_array,y_array,t);
u = interp1(t_array,u_array,t);
v = (24 + 24*a3*(b-t) + 12*a2*(b-t).^2 + 4*a1*(b-t).^3 + a0*(b-t).^4).*y...
    - ((b-t).^4).*u;
end