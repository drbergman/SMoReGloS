% This script will take the summary of the cohort data in summary.mat and
% add more summaries, such as mean across samples and standard deviation.

clearvars;
% addpath("~/Documents/MATLAB/myfunctions")

load("data/summary.mat","N","vals","cohort_size")

n_time_series = 1; % just the one time series of total cell counts
temp_avg = mean(N,ndims(N));
temp_std = std(N,[],ndims(N));
temp_avg = reshape(temp_avg,301,[]);
temp_std = reshape(temp_std,301,[]);
for i = size(temp_avg,2):-1:1
    D(i).A = temp_avg(:,i);
    D(i).S = temp_std(:,i);
end
D = reshape(D,[1,cohort_size]); % only one condition (e.g. not varying a dose)

par_names = ["progenitor division rate","stem cell division rate","tip migration rate","progenitor division limit"];
display_par_names = ["p_{div}","s_{div}","r_{mig}","p_{lim}"];
par_val_str = cellfun(@(c) arrayfun(@(i) string(num2str(c(i))),1:numel(c)),vals,'UniformOutput',false);

t = (0:300)*6/24; % time in hours

save("data/summary.mat","D","*par*","t","n_time_series","-append");

