a = [1:100000];
b = [1:100000];
res = [1:100000];
% Conventional Method
tic
for i = 1:100000
    res(i) = a(i)*b(i);
end
t1 = toc;
% Vectorization Method
res = a.*b;
t2= toc;

fprintf(" Conventional Method: %f\n", t1);
fprintf("Vectorization Method: %f\n", t2-t1);