function v = myfun_RK_25(t, t_array, y_array, t_i, a)
y = interp1(t_array,y_array,t);
v = ( -(t_i-t).*((t-a).^4) + 4*((t_i-t).^2).*((t-a).^3) - 2*((t_i-t).^3).*((t-a).^2) ).*y ;
end