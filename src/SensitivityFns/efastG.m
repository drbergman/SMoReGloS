function x = efastG(s,omega,phi)

% periodic functoin that gives uniform coverage on [0,1]
x = 0.5 + 0.318309886183791 * asin(sin(omega.*s+phi));
