function  [dx, ut, t_arr] = controlled_sys_sync_gen(t, x, t_array, R, R_dot, R_ddot, R_tdot, l0, l1, l2,b1,b2,b3,b4,E,P)
persistent u_array t_pts;

if(~isempty(x))
    dx = zeros(3,1);
%     b1 = 34.29; b2=0; b3=0.1490; b4=0.3341; P=28.220; E=0.2405;
    
    x1 = x(1); x2 = x(2); x3 = x(3);
    F_dot = x2; F_ddot = P -b2*F_dot - x3*b1*sin(x1);
    
    r = interp1(t_array, R, t); % interploted values
    r_dot = interp1(t_array, R_dot, t);
    r_ddot = interp1(t_array, R_ddot, t);
    r_tdot = interp1(t_array, R_tdot, t);
    
    v = r_tdot - l2*(F_ddot - r_ddot)...
        - l1*(F_dot - r_dot) - l0*(x1 - r);
    
    u = E - b4*x3 + b3*cos(x1)...
        - ((-b2*F_ddot - v)*sin(x1)-(P-b2*F_dot-F_ddot)*F_dot*cos(x1))/(b1*(sin(x1))^2);
    
    u_array = [u_array; u];
    t_pts = [t_pts; t];
    dx(1) = x2;
    dx(2) = -b1*x3*sin(x1) - b2*x2 + P;
    dx(3) = b3*cos(x1) - b4*x3 + E - u;
end

if nargout>1
    dx = -1;
    ut = u_array;
    t_arr = t_pts;
end
end