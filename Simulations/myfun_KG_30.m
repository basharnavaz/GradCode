function v = myfun_KG_30(t, t_array, y_array, b, a0, a1, a2)
y = interp1(t_array,y_array,t);
v = ( (6 + 6*a2.*(b-t) + 3*a1.*(b-t).^2 + a0.*(b-t).^3) ).*y ;
end