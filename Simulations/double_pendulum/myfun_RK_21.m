function v = myfun_RK_21(t, t_array, y_array, t_i, a)
y = interp1(t_array,y_array,t);
v = -( (1/6)*((t_i-t).^3).*((t-a).^4) ).*y ;
end
