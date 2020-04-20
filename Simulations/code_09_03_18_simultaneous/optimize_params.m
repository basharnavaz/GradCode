function [ true_p, diff ] = optimize_params(p,y_data,y_dot_data,y_ddot_data,y_tdot_data,U_data,phase,diff )

p0 = [34.29; 0; 0.1490; 0.3341; 28.220; 0.2405];
p_old = p;
opt_func = @(p)opt_func_SG_params(p, y_data, y_dot_data, y_ddot_data,...
            y_tdot_data, U_data, phase);
options = optimoptions('fmincon', 'MaxIterations', 10000,...
            'MaxFunEvals', 10000, 'UseParallel', true, 'Display', 'iter-detailed',...
            'Algorithm', 'sqp');
A=[]; b=[]; Aeq=[]; beq=[]; nonlcon = [];
lb = zeros(6,1);
ub=[100;100;100;100;100;100];
        
[p,fval] = fmincon(opt_func,p,A,b,Aeq,beq,lb,ub,nonlcon,options);
 
% update params
diff_new = norm(p-p0);
if(diff_new < diff)
    true_p = p;
    diff = diff_new; 
else
    true_p = p_old;
end

end

