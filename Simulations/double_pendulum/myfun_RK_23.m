function v = myfun_RK_23(t, t_array, y_array, t_i, a)
y = interp1(t_array,y_array,t);
v = ( -0.5*((t_i-t).^2).*((t-a).^4) + (2/3)*((t_i-t).^3).*((t-a).^3) ).*y ;
end
