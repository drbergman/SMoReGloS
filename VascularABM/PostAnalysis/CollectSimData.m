% This script will create a single .mat file of the cell counts in (all?)
% simulations. All simulations for now since I want all the simulations
% that I have run.

clearvars;

summary_file = "summary.mat";
vals = {[0.05,0.125,0.245], [0.01,0.05,0.1], [1,2,3],[8,12,16]};

n_conditions = 1; % number of conditions under which ABM was conducted (e.g. dosing values)
nsamps_per_parameter_vector = 6;
C = {[]}; % one condition, with no data to keep
% f = dir("../src/Data/Binary/*/*/fig_data/NumberofCells*t300.txt");
f = dir("../src/Data/Binary/*/*/fig_data/NumberofCells*t10.txt");
N = zeros(300,3,3,3,3,nsamps_per_parameter_vector); % [300 time points (after t=0), [4 pars at 3 vals each], 6 samples per par]
for i = 1:numel(f)
    I = findIndices(f(i).folder,vals);
    N(:,I{:}) = readmatrix(sprintf("%s/%s",f(i).folder,f(i).name));
end

N = cat(1,100+zeros([1,size(N,2:ndims(N))]),N);
cohort_size = size(N,2:(ndims(N)-1));

path_to_summary = sprintf("data/%s", summary_file);

[~,~,~] = mkdir("data");
save(path_to_summary,"N","vals","cohort_size","nsamps_per_parameter_vector","n_conditions","C");

function I = findIndices(F,V)

start_pats = {'AA_','AB_','BB_','CC_','Sample'};
d_to_start = [3,3,3,3,6];
end_pats = {'__AB','/BB','__CC','__Sample','/fig_data'};
I = cell(length(start_pats),1);
for i = 1:length(start_pats)
    i1 = strfind(F,start_pats{i});
    i2 = strfind(F,end_pats{i});
    if length(i1)~=1 || length(i2)~=1
        error("Wrong number of patterns found.")
    end
    str_val = F((i1+d_to_start(i)):(i2-1));
    str_val = regexprep(str_val,'_','\.');
    val = str2double(str_val);
    if strcmp(start_pats{i},'Sample')
        I{i} = val;
    else
        I{i} = find(V{i}==val);
    end
end
end