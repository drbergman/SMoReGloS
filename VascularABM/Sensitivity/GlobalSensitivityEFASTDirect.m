% This script will collect the HPC-run sims that will be used to compute
% global sensitivity.

clearvars;

%% identify the data
f = dir("data/sims/*eFAST*/Number*.txt");
load("data/EFASTSample.mat","points"); % points at which the MOAT samples were taken
load("../PostAnalysis/data/summary.mat","display_par_names")

%% what is the endpoint?
endpoints = ["final_size", "AUC", "time_to_half"];

for endpoint = endpoints
    t = 0:0.25:75;
    dt = t(2)-t(1);

    %% process points
    [nfacs,Ns,Nr,~] = size(points);


    %% process sim data
    N = zeros(numel(f),1); % just get the final tumor size
    for i = 1:numel(f)
        % Don't assume that dir returns the sims in the order I want
        I = strfind(f(i).folder,"_");
        num = str2double(f(i).folder(I(end)+1:end));
        assert(num==i);
        % now I know this one is in the proper spot in the order

        temp = readmatrix(sprintf("%s/%s",f(i).folder,f(i).name));
        if length(temp) < 300 || temp(300)==0
            error("This sim did not make it to the endpoint.")
        end
        
        switch endpoint
            case "final_size"
                N(i) = temp(end);
            case "AUC"
                N(i) = trapz(t,[100;temp]);
            case "time_to_half"
                temp = [100;temp];
                N_final = temp(end);
                if N_final <= 200
                    N(i) = 0;
                    continue;
                end
                t2 = find(temp>=0.5*N_final,1);
                t1 = t2-1;

                prop_dt = (0.5*N_final - temp(t1))/(temp(t2)-temp(t1));
                N(i) = t(t1) + prop_dt * dt;
        end
    end

    %% reshape and compute average
    N = reshape(N,[size(points,2:4),6]); % final dimension is for samples (I first loop over all parameter values, then over samples)
    A = mean(N,ndims(N));
    S = std(N,[],ndims(N));

    %% compute Fourier coefficients
    % copy details from MakeeFASTSamplePars.m
    omega_max = 8;
    M = 4;
    AA = A - mean(A,1);
    NQ = (Ns-1)/2;
    N0 = NQ+1;
    s = linspace(-pi,pi,Ns+1);
    s(end) = [];

    S1 = zeros(nfacs,1);
    ST = zeros(nfacs,1);

    for i=1:nfacs %loop through parameters
        for ri=1:Nr
            complement_total = 0;
            Y_VECP = AA((N0+1):end,ri,i)+AA((N0-1):-1:1,ri,i); % for convenience below and using cos(-theta) = cos(theta)
            Y_VECM = AA((N0+1):end,ri,i)-AA((N0-1):-1:1,ri,i); % for convenience below and using sin(-theta) = -sin(theta)
            for j=1:floor(omega_max/2)
                theta = j*s(N0+1:end);
                cos_theta = cos(theta);
                sin_theta = sin(theta);
                a_j(j) = (AA(N0,ri,i)+Y_VECP'*cos_theta')/Ns;
                b_j(j) = Y_VECM'*sin_theta'/Ns;
                complement_total = complement_total+a_j(j)^2+b_j(j)^2;
            end
            % Computation of V_{(ci)}.
            Vci(ri) = 2*complement_total;
            % Fourier coeff. at [P*omega_max, for P=1:M].
            first_order = 0;
            for j=omega_max:omega_max:omega_max*M
                % ANGLE = j*2*(1:NQ)*pi/Ns;
                theta = j*s(N0+1:end);
                cos_theta = cos(theta');
                sin_theta = sin(theta');
                a_j(j) = (AA(N0,ri,i)+Y_VECP'*cos_theta)/Ns;
                b_j(j) = Y_VECM'*sin_theta/Ns;
                first_order = first_order+a_j(j)^2+b_j(j)^2;
            end
            % Computation of V_i.
            Vi(ri) = 2*first_order;
            % Computation of the total variance
            % in the time domain.
            V(ri) = AA(:,ri,i)'*AA(:,ri,i)/Ns;
        end

        S1(i) = mean(Vi)/mean(V);
        ST(i,1,1) = 1-mean(Vci)/mean(V);
    end

    save("data/GlobalSensitivityEFASTDirect_" + endpoint,"M","Nr","Ns","S1","ST","display_par_names","omega_max");
end