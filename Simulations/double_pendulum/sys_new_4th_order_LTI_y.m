function dy = sys_new_4th_order_LTI_y(t, y)
dy = zeros(4,1);
dy(1) = y(2);
dy(2) = y(3);
dy(3) = y(4);
dy(4) = -156.96*y(3) - 1154.8332*y(1);
end