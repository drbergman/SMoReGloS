function out = computeSMAUC(p,~,model_type)

switch model_type
    case "exponential"
        out = 100 * (exp(75*p) - 1) / p;

    case "logistic"
        % r = p(1);
        % K = p(2);
        % out = y0*p(2)./((p(2)-y0).*exp(-p(1).*t_final)+y0);
        t = 0:.1:75;
        y = 100*p(2)./((p(2)-100).*exp(-p(1).*t)+100);
        out = trapz(t,y);

    case "von_bertalanffy"
        % compute log of endpoint since points in SM parameter space can
        % result in overflow errors
        p(2) = 1 - 1/p(2); % p(2) is nu coming in and theta going forward
        opts = odeset("NonNegative",1,"Events",@vbEvents);
        t = 0:.1:75;
        sol = ode45(@(t,x) odefn(x,p),[0 75],100,opts);

        if isempty(sol.ie)
            y = deval(sol,t);
            out = trapz(t,y);
            return;
        end

        %% otherwise deal with the events
        log_E = (1/(1-p(2))) * (log(p(1))-log(p(3))); % set it to the (log of) the equilibrium value
        if log_E < 0
            out = 7500/(log_E-log(100))*(exp(log_E)/100-1);
            out = log(out);
        else
            out = log(7500/(log_E-log(100))) + log_E - log(100); % close to the AUC
        end
    otherwise
        error("%s is an unspecified SM model.\n",model_type);

end