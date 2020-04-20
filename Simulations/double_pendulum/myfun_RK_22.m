function v = myfun_RK_22(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
v = ( (1/6)*((t_i-t).^3).*((b-t).^4) ).*y ;
end