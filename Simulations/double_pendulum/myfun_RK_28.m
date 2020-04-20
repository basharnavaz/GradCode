function v = myfun_RK_28(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
v = ( (b-t).^4 + 12*(t_i-t).*((b-t).^3) + 18*((t_i-t).^2).*((b-t).^2) + 4*((t_i-t).^3).*(b-t) ).*y ;
end