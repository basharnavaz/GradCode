function [v] = myfun_KG_26(t, t_array, y_array, t_i, b,a0,a1,a2)
y = interp1(t_array,y_array,t);
v = ( 9*((b-t).^2 + a2*(b-t).^3) + (t_i-t).*(18*(b-t) + 6*a2.*(b-t).^2 + a1.*((b-t).^3)) + 0.5*((t_i-t).^2).*(6 + 6*a2.*(b-t) + 3*a1.*((b-t).^2) + a0.*(b-t).^3) ).*y ;
end

