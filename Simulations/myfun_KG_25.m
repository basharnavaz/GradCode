function [v] = myfun_KG_25(t, t_array, y_array, t_i, a, a0, a1, a2)
y = interp1(t_array,y_array,t);
v = ( 9*((t-a).^2 - a2.*(t-a).^3) +(t_i-t).*(-18*(t-a) + 6*a2.*(t-a).^2 -a1.*(t-a).^3) + 0.5*((t_i-t).^2).*(6 - 6*a2.*(t-a) + 3*a1.*((t-a).^2)-a0.*(t-a).^3) ).*y ;
end

