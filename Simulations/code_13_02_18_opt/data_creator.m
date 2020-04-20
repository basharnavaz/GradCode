l1 = y_tdot_true(1:1199);
l2 = y_tdot_data(1:1199)';
diff=l1-l2;
good_estimates = diff>0 & diff<1;
sum(good_estimates);
good_estimate_indices = find(good_estimates~=0);
matched_ests = zeros(length(good_estimate_indices),8);
matched_ests(:,1) = y_true(good_estimate_indices)';
matched_ests(:,2) = y_data(good_estimate_indices);
matched_ests(:,3) = y_dot_true(good_estimate_indices)';
matched_ests(:,4) = y_dot_data(good_estimate_indices);
matched_ests(:,5) = y_ddot_true(good_estimate_indices)';
matched_ests(:,6) = y_ddot_data(good_estimate_indices);
matched_ests(:,7) = y_tdot_true(good_estimate_indices)';
matched_ests(:,8) = y_tdot_data(good_estimate_indices);

save('estimation_win50_control_1.mat')