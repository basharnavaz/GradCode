function v = myfunc_NKM_2(t, t_array, y_array, t_i, b,a0,a1,a2,a3)
y = interp1(t_array,y_array,t);
%for K_b
v=((72.*(b-t).^2 + 12.*a2.*(b-t).^3 + a2.*(b-t).^4) + (t-t_i).*(96.*(b-t)+ 36.*a3.*(b-t).^2 + 8.*a2.*(b-t).^3 + a1.*(b-t).^4) + 0.5.*(t-t_i).^2.*(24 + 24.*a3.*(b-t) + 12.*a2.*(b-t) + 12.*a2.*(b-t).^2 + 4.*a2.*(b-t).^3 + a0.*(b-t).^4)).*y;
end