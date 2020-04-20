function [obj_func_value] = opt_func_sync_gen(X, y, y_dot, y_ddot, y_tdot, u)

    b1 = X(1);
    b2 = X(2);
    b3 = X(3);
    b4 = X(4);
    P = X(5); 
    E = X(6); 

    sin_y = sin(y);
    cos_y = cos(y);
    vector_func = (E*b1*(sin_y).^2 - b4*P.*(sin_y) + b2*b4.*y_dot.*sin_y...
                + b4*y_ddot.*sin_y + b1*b3.*(sin_y.^2).*cos_y...
                + b2*(sin_y.*y_ddot - (y_dot.^2).*cos_y)...
                + P.*y_dot.*cos_y - b1.*sin_y.^2.*u ...
                + y_tdot.*sin_y - y_dot.*y_ddot.*cos_y).^2;

    obj_func_value = sum(vector_func);
%     obj_func_value = -obj_func_value;

end

