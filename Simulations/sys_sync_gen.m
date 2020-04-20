function  dx = sys_sync_gen(t, x, u_array, t_array)
   dx = zeros(3,1);
   x1 = x(1); x2=x(2); x3=x(3);
%  using true parameters
   b1 = 34.29; b2=0; b3=0.1490; b4=0.3341; P=28.220; E=0.2405;
   u = interp1(t_array,u_array,t);
   dx(1) = x2;
   dx(2) = -b1*x3*sin(x1) - b2*x2 + P;
   dx(3) = b3*cos(x1) - b4*x3 + E - u;
end