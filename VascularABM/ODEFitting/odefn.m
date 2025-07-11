function dx = odefn(x,p)

% the ODE function defining the SM

dx = p(1)*x.^(p(2)) - p(3)*x;
