function out = computeSMEndpoint(p,~,model_type)

% y0 = 100;
switch model_type
    case "exponential"
        out = 100 * exp(75*p);

    case "logistic"
        % r = p(1);
        % K = p(2);
        % out = y0*p(2)./((p(2)-y0).*exp(-p(1).*t_final)+y0);
        out = 100*p(2)./((p(2)-100).*exp(-p(1).*75)+100);

    case "von_bertalanffy"
        % compute log of endpoint since points in SM parameter space can
        % result in overflow errors
        p(2) = 1 - 1/p(2); % p(2) is nu coming in and theta going forward
        opts = odeset("NonNegative",1,"Events",@vbEvents);
        sol = ode45(@(t,x) odefn(x,p),[0,75],100,opts);

        if isempty(sol.ie)
            out = log(sol.y(end));
            return;
        end

        %% otherwise deal with the events
        out = (1/(1-p(2))) * (log(p(1))-log(p(3))); % set it to the (log of) the equilibrium value

    otherwise
        error("%s is an unspecified SM model.\n",model_type);

end