function [obj_func_value] = opt_func_SG_params(p, y, y_dot, y_ddot, y_tdot, u, flag)

    b1 = p(1);
    b2 = p(2);
    b3 = p(3);
    b4 = p(4);
    P = p(5); 
    E = p(6); 
    
%     if(length(y)>50)
%         y = y(end-50 : end);
%         y_dot = y_dot(end-50 : end);
%         y_ddot = y_ddot(end-50 : end);
%         y_tdot = y_tdot(end-50 : end);
%         u = u(end-50 : end);
%     end

    sin_y = sin(y);
    cos_y = cos(y);
    vector_func = (E*b1*(sin_y).^2 - b4*P.*(sin_y) + b2*b4.*y_dot.*sin_y...
                + b4*y_ddot.*sin_y + b1*b3.*(sin_y.^2).*cos_y...
                + b2*(sin_y.*y_ddot - (y_dot.^2).*cos_y)...
                + P.*y_dot.*cos_y - b1.*sin_y.^2.*u ...
                + y_tdot.*sin_y - y_dot.*y_ddot.*cos_y).^2;
%     if flag == 1
        obj_func_value = sum(vector_func) + norm(p)^2;
%     else
%         obj_func_value = sum(vector_func);
%     end

end

