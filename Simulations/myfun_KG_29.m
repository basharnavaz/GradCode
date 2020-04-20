function v = myfun_KG_29(t, t_array, y_array, a, a0, a1, a2)
y = interp1(t_array,y_array,t);
v = ( (6 - 6*a2.*(t-a) + 3*a1.*(t-a).^2 - a0.*(t-a).^3) ).*y ;
end