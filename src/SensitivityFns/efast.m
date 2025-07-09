% function [S1,ST,order] = efast(f,n,Nr,omega_max,M,Ns,options)
function [S1,ST] = efast(f,n,Nr,omega_max,M,Ns)

% f = studied function
% n = number of factors
% Nr = number of resamplings for each factor

arguments
    f
    n double {mustBeInteger}
    Nr double {mustBeInteger} = 1
    omega_max double {mustBeInteger} = 8
    M double {mustBeInteger} = 4;
    Ns double {mustBeInteger} = 2*M*omega_max + 1
    % options.sort_output_by_st logical = true
end

% equally-spaced points on (-pi,pi) but only include left endpoint
s = linspace(-pi,pi,Ns+1);
s(end) = [];


%% compute values along s for each factor and each resample
y = zeros(Ns,n,Nr);
for fi = 1:n % factor index
    omega = ones(1,n);
    omega(fi) = omega_max;

    % set the complementary frequencies to the values between 1 and floor(omega/2)
    omega(setdiff(1:n,fi)) = mod(0:n-2,floor(omega_max/2)) + 1;

    for ri = 1:Nr
        phi = 2*pi*rand(1,n); % random phase shifts
        x = efastG(s,omega',phi');

        parfor si = 1:Ns
            y(si,fi,ri) = f(x(:,si));
        end

    end
end

%% compute Fourier coefficients
Y = y - mean(y,1);
NQ = (Ns-1)/2;
N0 = NQ+1;
complement_var = zeros(Nr,n);
first_order_var = zeros(Nr,n);
total_var = zeros(Nr,n);
for fi=1:n
    for ri=1:Nr
        complement_total = 0;
        Y_VECP = Y((N0+1):end,fi,ri)+Y((N0-1):-1:1,fi,ri); % for convenience below and using cos(-theta) = cos(theta)
        Y_VECM = Y((N0+1):end,fi,ri)-Y((N0-1):-1:1,fi,ri); % for convenience below and using sin(-theta) = -sin(theta)
        a_j = zeros(1,floor(omega_max/2));
        b_j = zeros(1,floor(omega_max/2));
        for j=1:floor(omega_max/2)
            theta = j*s(N0+1:end);
            cos_theta = cos(theta);
            sin_theta = sin(theta);
            a_j(j) = (Y(N0,fi,ri)+Y_VECP'*cos_theta')/Ns;
            b_j(j) = Y_VECM'*sin_theta'/Ns;
            complement_total = complement_total+a_j(j)^2+b_j(j)^2;
        end
        complement_var(ri,fi) = 2*complement_total;
        % Fourier coeff. at [P*omega_max, for P=1:M].
        first_order = 0;
        a_j = zeros(1,length(omega_max:omega_max:omega_max*M));
        b_j = zeros(1,length(omega_max:omega_max:omega_max*M));
        for j=omega_max:omega_max:omega_max*M
            % ANGLE = j*2*(1:NQ)*pi/Ns;
            theta = j*s(N0+1:end);
            cos_theta = cos(theta');
            sin_theta = sin(theta');
            a_j(j) = (Y(N0,fi,ri)+Y_VECP'*cos_theta)/Ns;
            b_j(j) = Y_VECM'*sin_theta/Ns;
            first_order = first_order+a_j(j)^2+b_j(j)^2;
        end
        first_order_var(ri,fi) = 2*first_order;
        total_var(ri,fi) = Y(:,fi,ri)'*Y(:,fi,ri)/Ns;
    end
end

S1 = mean(first_order_var,1)./mean(total_var,1);
ST = 1-mean(complement_var,1)./mean(total_var,1);

% if options.sort_output_by_st
%     [ST,order] = sort(ST,"descend");
%     S1 = S1(order);
% else
%     order = 1:length(S1);
% end


end
