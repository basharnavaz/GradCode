function v = myfunc_NKM_1(t, t_array, y_array, t_i, a,a0,a1,a2,a3)
y = interp1(t_array,y_array,t);
%for K_f
v = ((72.*(t-a).^2 - 12.*a3.*(t-a).^3 + a2.*(t-a).^4)+(t-t_i).*(-96.*(t-a) + 36.*a3.*(t-a).^2 - 8.*a2.*(t-a).^3 + 4.*a1.*(t-a).^4)+0.5.*(t-t_i).^2.*(24 - 24.*a3.*(t-a) + 12.*a2.*(t-a).^2 - 4.*a1.*(t-a).^3 + a0.*(t-a).^4)).*y;
end