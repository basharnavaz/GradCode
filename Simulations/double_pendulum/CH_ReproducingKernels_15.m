% KG 

clc
clear all

fs = 1000; % sampling rate

%system information%
a = 0;
b = 2;
% x0 = [(0.0523599); (0.0872665); 0; 0];
% t_array = a:(1/fs):b;
% [t_array, x_array] = ode45(@sys_DoublePendulum,t_array,x0);
% x_array = x_array';
% y = [0 1 0 0]*x_array;
% y = y';
x0 = [(0.0523599); (0.0872655); 0; 0];
t_array = a:(1/fs):b;
[t_array, x_array] = ode45(@sys_new_4th_order_LTI_y,t_array,x0);
x_array = x_array';
y = [0 1 0 0]*x_array;
y = y';
plot(t_array,y)

%addition of noise%
y_noisy = awgn(y,37,'measured');%measured signal with error%
y_noisy = y_noisy';
plot(t_array,y_noisy)
% y_n = y_noisy;
y_n = y;

V5 = 0;
V6 = 0;
V7 = 0;
V8 = 0;
V9 = 0;
V10 = 0;
V11 = 0;
V12 = 0;
V13 = 0;
V14 = 0;
V15 = 0;
V16 = 0;
V17 = 0;
V18 = 0;
V19 = 0;
V20 = 0;
V21 = 0;
V22 = 0;
V23 = 0;
V24 = 0;

max_count = 100000; % Max interval count for quadgk

Num = 0; 

%Parameter estimation using K_DS 

for i = 1:1:length(t_array)
    
    t_i = t_array(i);

    Num = Num+1;
   
    W0 = quadgk(@(t,v)myfun_RK_21(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_22(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    W1 = quadgk(@(t,v)myfun_RK_23(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_24(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    W2 = quadgk(@(t,v)myfun_RK_25(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_26(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    W3 = quadgk(@(t,v)myfun_RK_27(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_28(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    V0 = quadgk(@(t,v)myfun_RK_29(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_30(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    V1 = quadgk(@(t,v)myfun_RK_31(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_32(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    V2 = quadgk(@(t,v)myfun_RK_33(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_34(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
    V3 = quadgk(@(t,v)myfun_RK_35(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_36(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    
%     V4 = -((t_i-a).^4+(b-t_i).^4)*y(i) + quadgk(@(t,v)myfun_RK_37(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_38(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;
    V4 = -((t_i-a).^4+(b-t_i).^4)*y_n(i) + quadgk(@(t,v)myfun_RK_37(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_38(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

    V5 = V5 + V0*W0;
    V6 = V6 + V1*W0;
    V7 = V7 + V2*W0;
    V8 = V8 + V3*W0;
    V9 = V9 + V4*W0;
    V10 = V10 + V0*W1;
    V11 = V11 + V1*W1;
    V12 = V12 + V2*W1;
    V13 = V13 + V3*W1;
    V14 = V14 + V4*W1;
    V15 = V15 + V0*W2;
    V16 = V16 + V1*W2;
    V17 = V17 + V2*W2;
    V18 = V18 + V3*W2;
    V19 = V19 + V4*W2;
    V20 = V20 + V0*W3;
    V21 = V21 + V1*W3;
    V22 = V22 + V2*W3;
    V23 = V23 + V3*W3;
    V24 = V24 + V4*W3;
    
end

G1 = [V5 V6 V7 V8; V10 V11 V12 V13; V15 V16 V17 V18; V20 V21 V22 V23];
G2 = [-V9; -V14; -V19; -V24];
% G1_inv = inv(G1);
aa = G1\G2;
a0 = aa(1)
a1 = aa(2)
a2 = aa(3)
a3 = aa(4)
plot(t_array,y_noisy)
