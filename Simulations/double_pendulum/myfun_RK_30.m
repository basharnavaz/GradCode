function v = myfun_RK_30(t, t_array, y_array, t_i, b)
y = interp1(t_array,y_array,t);
% v = ( 4 + 16*(b.^3) + 48*b.*t - 120*(b.^2).*t - 48*(t.^2) + 192*b.*(t.^2) - 88*(t.^3) - 48*b.*t_i + 72*(b.^2).*t_i + 48*t.*t_i - 144*b.*t.*t_i + 72*(t.^2).*t_i ).*y;
% v = ( (1/6).*(b.^4) - (2/3).*(b.^3).*t + (b.^2).*(t.^2) - (2/3)*b.*(t.^3) + (1/6)*(t.^4) ).*y;
v = ( -(1/6).*(b.^4).*(t.^3) + (2/3).*(b.^3).*(t.^4) - (b.^2).*(t.^5) + (2/3)*b.*(t.^6) - (1/6)*(t.^7) + (1/2)*(b.^4).*(t.^2).*t_i - 2*(b.^3).*(t.^3).*t_i + 3*(b.^2).*(t.^4).*t_i - 2*b.*(t.^5).*t_i + (1/2)*(t.^6).*t_i - (1/2)*(b.^4).*t.*(t_i.^2) +  2*(b.^3).*(t.^2).*(t_i.^2) - 3*(b.^2).*(t.^3).*(t_i^2) + 2*b.*(t.^4).*(t_i.^2) - (1/2)*(t.^5).*(t_i.^2) + (1/6)*(b.^4).*(t_i.^3) - (2/3)*(b.^3).*t.*(t_i.^3) + (b.^2).*(t.^2).*(t_i.^3) - (2/3)*b.*(t.^3).*(t_i.^3) + (1/6)*(t.^4).*(t_i.^3) ).*y;
end