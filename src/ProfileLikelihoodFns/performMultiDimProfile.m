function MDProfile = performMultiDimProfile(files,objfn_constants,par_vals,input_opts)

opts = defaultPerformMultiDimProfileOptions;
if nargin >= 4 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

sz = reshape(cellfun(@numel,par_vals),1,[]);
N = prod(sz);
npars = length(sz);

load(files.data_file,"t","D","C","cohort_size");
m = size(D,1);
D = reshape(D,m,[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used


field_sz = size(D(1).A);
nt = field_sz(1);
n_time_series = field_sz(2);
n = prod(field_sz);
A = arrayify(D,"A",true);
A = reshape(A,[field_sz,m,n_abm_vecs]);
if opts.assume_independent_time_series
    S = arrayify(D,"S",true);
    S = reshape(S,[field_sz,m,n_abm_vecs]);
else
    dets = arrayify(D,"D",true);
    dets = reshape(dets,[nt,m,n_abm_vecs]);
    Q = arrayify(D,"I",true);
    Q = reshape(Q,[n_time_series,n_time_series,nt,m,n_abm_vecs]);
end

fn = objfn_constants.fn;
fn_opts = objfn_constants.fn_opts;

P = allCombos(par_vals{:},'matlab');

if isfield(files,"previous_profile_file")
    error("Not yet sure how to reload and start from there.")
    load(files.previous_profile_file,"MDProfile") %#ok<UNRCH>
    MDProfile = reshape(MDProfile,n_abm_vecs,N);
else
    MDProfile = zeros(n_abm_vecs,N);
end

dim_constant = -0.5*n_time_series*log(2*pi());

if opts.force_serial
    if opts.assume_independent_time_series
        for i = 1:N
            p = P(i,:)';
            sim_data = zeros([field_sz,m]);
            for j = 1:m
                sim_data(:,:,j) = fn(p,t,C{j},fn_opts);
            end
            MDProfile(:,i) = LLIndependent(dim_constant,S,A,sim_data);
        end
    else
        sim_data = zeros([field_sz,m,N]);
        for i = 1:N
            p = P(i,:)';
            for j = 1:m
                sim_data(:,:,j,i) = fn(p,t,C{j},fn_opts);
            end

        end
        MDProfile(:,i) = MDProfile(:,i) - 0.5*reshape(sum(log(dets),1:2),[],1);
        for ti = 1:nt
            for k = 1:m
                for j = 1:n_abm_vecs
                    temp = reshape(sim_data(ti,:,k,:),n_time_series,N);
                    MDProfile(j,:) = MDProfile(j,:) - 0.5 * sum(temp .* (Q(:,:,ti,k,j) * temp));
                end
            end
        end
    end
else
    if isempty(gcp('nocreate'))
        ppool = parpool("Threads"); %#ok<NASGU>
    else
        ppool = gcp; %#ok<NASGU>
    end
    if opts.assume_independent_time_series
        C = parallel.pool.Constant(C);
        parfor i = 1:N
            p = P(i,:)';
            sim_data = zeros([field_sz,m]);
            for j = 1:m
                sim_data(:,:,j) = fn(p,t,C.Value{j},fn_opts); %#ok<PFBNS> (suppresses warning about fn bearing an overhead cost
            end
            MDProfile(:,i) = LLIndependent(dim_constant,S,A,sim_data);
        end
    else
        error("Haven't set this up yet")
    end
end

MDProfile = reshape(MDProfile,[cohort_size,sz]);

end

function default_options = defaultPerformMultiDimProfileOptions

default_options.force_serial = true; % whether to allow for parallelization
default_options.assume_independent_time_series = true; % assume that the time series are independent

end

function out = LLIndependent(dim_constant,S,A,sim_data)
out = dim_constant - 0.5*(sum(S.^2,1:3) + sum(((sim_data - A)./S).^2,1:3,"omitnan"));
end
