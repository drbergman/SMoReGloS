function p = basePars(model_type)


switch model_type
    case "exponential"
        p = 0.05; % lambda

    case "logistic"
        % This load the base logistic parameters for x' = r * x * (1-x/K)
        p = zeros(2,1);
        p(1) = 1e-2; % r
        p(2) = 1e4; % K

    case "von_bertalanffy"
        % this loads up what I will consider the base parameter values for the SM
        % This SM is given by x' = alpha * x^theta - beta * x

        p = zeros(3,1);

        p(1) = 1; % alpha
        p(2) = 2; % nu (theta = 1 - 1/nu) (this allows nu to be on the open interval (0,inf) which hopefully helps with optimization)
        p(3) = 1; % beta

    otherwise
        error("%s is an unspecified SM model.\n",model_type);
end