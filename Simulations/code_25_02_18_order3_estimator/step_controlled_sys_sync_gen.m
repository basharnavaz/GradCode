function  [dx] = step_controlled_sys_sync_gen(t, x, u)

dx = zeros(3,1);

x1 = x(1); x2 = x(2); x3 = x(3);
b1 = 34.29; b2=0; b3=0.1490; b4=0.3341; P=28.220; E=0.2405;
dx(1) = x2;
dx(2) = -b1*x3*sin(x1) - b2*x2 + P;
dx(3) = b3*cos(x1) - b4*x3 + E - u;

end