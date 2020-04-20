X0 = ones(6,1);
% p=[11.4548 0.9231 15.316 -0.3081 28.5004 11.4551];
% X0 = [34 0.1 0.1 0.1 0.1 25];
X0 = [50 50 50 50 50 50];
X0 = [-1 2 3 4 6 7]
% X0 = [1000 1000 1000 1000 1000 1000];
C0 = 340; 
% X0 = [34.28; 0; 0.1480; 0.3340; 28.210; 0.2400];
% X0=[34.29; 0; 0.1490; 0.3341; 28.220; 0.2405];
%  TRUE PARAMS: b1 = 34.29; b2=0; b3=0.1490; b4=0.3341; P=28.220; E=0.2405;
y=matched_ests(:,2);
y_dot=matched_ests(:,3);
y_ddot=matched_ests(:,4);
y_tdot=matched_ests(:,5);
u=matched_ests(:,6);

opt_f = @(X)opt_func_sync_gen(X0, y, y_dot, y_ddot, y_tdot, u);
options = optimoptions('fmincon', 'MaxIterations', 10000,...
            'MaxFunEvals', 10000, 'UseParallel', true, 'Algorithm','sqp');
A=[];
b=[];
Aeq=[];
beq=[];
lb = zeros(6,1);
ub=[100;100;100;100;100;100];
nonlcon = [];
[Y,fval] = fmincon(opt_f,X0,A,b,Aeq,beq,lb,ub,nonlcon,options);

b1 = Y(1);
b2 = Y(2);
b3 = Y(3);
b4 = Y(4);
P = Y(5);
E = Y(6);
fprintf('%.3f %.3f %.3f %.3f %.3f %.3f %.3f',b1,b2,b3,b4,P,E);