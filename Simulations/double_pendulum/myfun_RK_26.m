function v = myfun_RK_26(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
v = ( (t_i-t).*((b-t).^4) + 4*((t_i-t).^2).*((b-t).^3) + 2*((t_i-t).^3).*((b-t).^2) ).*y ;
end