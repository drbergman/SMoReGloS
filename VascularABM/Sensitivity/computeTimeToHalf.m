function out = computeTimeToHalf(p,model_type)

switch model_type
    case "exponential"
        out = 75 - log(2)/p;

    case "logistic"
        % r = p(1);
        % K = p(2);
        out = -log((2*(p(2)-100)*exp(-75*p(1)) + 100)/(p(2)-100))/p(1); % time to half the max value

    case "von_bertalanffy"
        error("Please don't make me do this.")
   
    otherwise
        error("%s is an unspecified SM model.\n",model_type);

end