function v = myfun_RK_34(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
% v = ( 2*(b.^2) - 4*b.*t + 4*(b.^3).*t - (b.^4).*t + 2*(t.^2) - 12*(b.^2).*(t.^2) + 4*(b.^3).*(t.^2) + 12*b.*(t.^3) - 6*(b.^2).*(t.^3) - 4*(t.^4) + 4*b.*(t.^4) - (t.^5) - 4*(b.^3).*t_i + (b.^4).*t_i + 12*(b.^2).*t.*t_i - 4*(b.^3).*t.*t_i - 12*b.*(t.^2).*t_i + 6*(b.^2).*(t.^2).*t_i + 4*(t.^3).*t_i - 4*b.*(t.^3).*t_i + (t.^4).*t_i ).*y;
v = ( -(b.^4).*t + 8*(b.^3).*(t.^2) - 20*(b.^2).*(t.^3) + 20*b.*(t.^4) - 7*(t.^5) + (b.^4).*t_i - 12*(b.^3).*t.*t_i + 36*(b.^2).*(t.^2).*t_i - 40*b.*(t.^3).*t_i + 15*(t.^4).*t_i + 4*(b.^3).*(t_i.^2) - 18*(b.^2).*t.*(t_i.^2) + 24*b.*(t.^2).*(t_i.^2) - 10*(t.^3).*(t_i.^2) + 2*(b.^2).*(t_i.^3) - 4*b.*t.*(t_i.^3) + 2*(t.^2).*(t_i.^3) ).*y;
end