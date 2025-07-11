function out = computeTimeSeries(p,tt,~,fn_opts,~)

% y0 = 100;
switch fn_opts.model_type
    case "exponential"
        % lambda = p(1);
        % y' = lambda*y
        out = 100 * exp(p * tt');

    case "logistic"
        % r = p(1);
        % K = p(2);
        % out = y0*p(2)./((p(2)-y0).*exp(-p(1).*tt')+y0);
        out = 100*p(2)./((p(2)-100).*exp(-p(1).*tt')+100);
        out(tt==0) = 100; % just make sure the initial value is exactly 100

    case "von_bertalanffy"

        if isfield(fn_opts,"enforce_inequality") && fn_opts.enforce_inequality && p(3)>p(1) % beta should not be bigger than alpha (or else the solution will immediately start to decrease for any x>1 since theta<1)
            out = 100*exp((p(1)-p(3))*tt'); % just approximate this as exponential decay to avoid issues of the solution going negative and giving complex answers that even NonNegativity constraints cannot fix
            return;
        end
        p(2) = 1 - 1/p(2); % p(2) is nu coming in and theta going forward
        opts = odeset("NonNegative",1);
        sol = ode45(@(t,x) odefn(x,p),[0 tt(end)],100,opts);
        out = deval(sol,tt)';

        if any(isnan(out))
            out(isnan(out)) = -1e6; % something really different from whatever data I'm comparing it to so that it gets heavily penalized for getting NaN
        end

        if length(out)<length(tt) % for some parameter values, the solution blows up in finite time and this penalized those
            out((length(out)+1):length(tt)) = -1e6;
        end

        % out = min(out,1e10); % for plenty of chosen parameters (alpha >> beta for example) the solutions can grow large, but finite, s.t. the error they produce is Inf

    otherwise
        error("%s is an unspecified SM model.\n",fn_opts.model_type)
end
