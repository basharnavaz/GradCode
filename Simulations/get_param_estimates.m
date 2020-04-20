function [ a0, a1, a2 ] = get_param_estimates( t_array, a, b, y_n, max_count )
    V5 = 0; V6 = 0; V7 = 0; V8 = 0; V9 = 0; V10 = 0;
    V11 = 0; V12 = 0; V13 = 0; V14 = 0; V15 = 0; V16 = 0;


    for i = 1:100:length(t_array)

        t_i = t_array(i);

        W1 = quadgk(@(t,v)myfun_RK_1(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_2(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        W2 = quadgk(@(t,v)myfun_RK_3(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_4(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        W3 = quadgk(@(t,v)myfun_RK_5(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_6(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        V1 = quadgk(@(t,v)myfun_RK_7(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_8(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        V2 = quadgk(@(t,v)myfun_RK_9(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_10(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        V3 = quadgk(@(t,v)myfun_RK_11(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_12(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;

        V4 = ((t_i-a).^3+(b-t_i).^3)*y_n(i) + quadgk(@(t,v)myfun_RK_13(t,t_array,y_n,t_i, a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)myfun_RK_14(t,t_array,y_n,t_i, b),t_i,b,'MaxIntervalCount',max_count) ;


        V5 = V5 + V1*W1;
        V6 = V6 + V2*W1;
        V7 = V7 + V3*W1;
        V8 = V8 + V4*W1;
        V9 = V9 + V1*W2;
        V10 = V10 + V2*W2;
        V11 = V11 + V3*W2;
        V12 = V12 + V4*W2;
        V13 = V13 + V1*W3;
        V14 = V14 + V2*W3;
        V15 = V15 + V3*W3;
        V16 = V16 + V4*W3;

    end

    G1 = [V5 V6 V7; V9 V10 V11; V13 V14 V15];
    G2 = [-V8; -V12; -V16];
    % G1_inv = inv(G1);
    % aa = G1_inv*G2;
    aa = G1\G2;
    a0 = (-1)*aa(1); % multiplying by (-1) for a derivation sign mismatch
    a1 = (-1)*aa(2);
    a2 = (-1)*aa(3);

end

