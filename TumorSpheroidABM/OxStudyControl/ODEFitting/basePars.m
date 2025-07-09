function [p,lb,ub] = basePars()

p = zeros(3,1);
p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K

lb = [0;0;0];
ub = [10;20;1e4];
