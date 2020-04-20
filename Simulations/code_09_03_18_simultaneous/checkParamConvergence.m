function [ flag ] = checkParamConvergence(p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p0 = [34.29; 0;0.1490;0.3341;28.220; 0.2405];
diff = norm(p0-p);
if(diff<0.5)
    flag = 1;
else
    flag = 0;
end
end

