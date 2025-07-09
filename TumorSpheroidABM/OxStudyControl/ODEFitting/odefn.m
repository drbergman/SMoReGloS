function dx = odefn(~,x,p,~)

% the ODE function defining the SM

% dx = [-1;1] * p.lambda * x(1) + p.alpha * [1+p.P(x);-1] * x(2) - p.delta * x;
% dx = [-1;1] * p(1) * x(1) + p(2) * [2 - sum(x)/p(3);-1] * x(2) - p(4) * x;
dx = [-1;1] * p(1) * x(1) + p(2) * [2 - sum(x)/p(3);-1] * x(2); % for control, we got rid of the death rate
