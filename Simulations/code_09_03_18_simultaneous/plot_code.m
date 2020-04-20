% SIM 2

figure
subplot(2,1,1)
plot(t_array,y_real)
hold on
plot(t_array,y_est)
grid on
legend({'$\mathbf{y_{true}}$','$\mathbf{y_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
subplot(2,1,2)
plot(t_array,y_real)
hold on
plot(t_array,y_est)
grid on
legend({'$\mathbf{y_{true}}$','$\mathbf{y_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
%%
figure
subplot(2,1,1)
plot(t_array,y_dot_real)
hold on
plot(t_array,y_dot_est)
grid on
legend({'$\mathbf{y^{(1)}_{true}}$','$\mathbf{y^{(1)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
subplot(2,1,2)
plot(t_array,y_dot_real)
hold on
plot(t_array,y_dot_est)
grid on
legend({'$\mathbf{y^{(1)}_{true}}$','$\mathbf{y^{(1)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');

%%
figure
subplot(2,1,1)
plot(t_array,y_ddot_real)
hold on
plot(t_array,y_ddot_est)
grid on
legend({'$\mathbf{y^{(2)}_{true}}$','$\mathbf{y^{(2)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
subplot(2,1,2)
plot(t_array,y_ddot_real)
hold on
plot(t_array,y_ddot_est)
grid on
legend({'$\mathbf{y^{(2)}_{true}}$','$\mathbf{y^{(2)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');

%%

figure
subplot(2,1,1)
plot(t_array,y_tdot_real)
hold on
plot(t_array,y_tdot_est)
grid on
legend({'$\mathbf{y^{(3)}_{true}}$','$\mathbf{y^{(3)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');
subplot(2,1,2)
plot(t_array,y_tdot_real)
hold on
plot(t_array,y_tdot_est)
grid on
legend({'$\mathbf{y^{(3)}_{true}}$','$\mathbf{y^{(3)}_{estimated}}$'},'Interpreter','latex')
yL = get(gca,'YLim');
line([t_switch t_switch],yL,'Color','g');

%% SIM 1
y_est = interp1(x_data, y_data, t_array, 'spline');
plot(t_array,y_est)
hold on
plot(t_array,ref_traj)
legend({'Controlled output $\mathbf{y(t)}$', 'Reference Trajectory $\mathbf{F^{*}(t)}$'}, 'Interpreter','latex');
yL = get(gca,'YLim');
line([0.6085 0.6085],yL,'Color','g');
grid on
%
uu = interp1(x_data, U_data, t_array, 'spline');
plot(t_array,uu)
legend({'Control signal $\mathbf{u^{*}(t)}$'}, 'Interpreter', 'latex');
yL = get(gca,'YLim');
line([0.6085 0.6085],yL,'Color','g');
grid on
