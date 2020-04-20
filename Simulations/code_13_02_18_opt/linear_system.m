A = zeros(6,8);
b = zeros(6,1);
for i=1:8
    F = matched_ests(i,2);
    F_dot = matched_ests(i,3);
    F_ddot = matched_ests(i,4);
    F_tdot = matched_ests(i,5);
    
    t1 = sin(F)^2;
    t2 = -sin(F);
    t3 = F_dot*sin(F);
    t4 = F_ddot*sin(F);
    t5 = (sin(F)^2)*cos(F);
    t6 = sin(F)*F_ddot-(F_dot^2)*cos(F);
    t7 = F_dot*cos(F);
    t8 = -(sin(F)^2);
    a=[t1 t2 t3 t4 t5 t6 t7 t8];
    A(i,:) = a;
    b(i) = -F_tdot*sin(F) + F_dot*F_ddot*cos(F);
end