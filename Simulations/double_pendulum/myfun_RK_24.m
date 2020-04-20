function v = myfun_RK_24(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
v = ( 0.5*((t_i-t).^2).*((b-t).^4) + (2/3)*((t_i-t).^3).*((b-t).^3) ).*y ;
end