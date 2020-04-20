function v = myfun_RK_27(t, t_array, y_array, t_i, a)
y = interp1(t_array,y_array,t);
v = ( -((t-a).^4) + 12*(t_i-t).*((t-a).^3) - 18*((t_i-t).^2).*((t-a).^2) + 4*((t_i-t).^3).*(t-a) ).*y ;
end