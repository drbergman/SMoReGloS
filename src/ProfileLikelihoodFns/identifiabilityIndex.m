function indices = identifiabilityIndex(files)

load(files.profiles,"profiles")

profiles = reshape(profiles,size(profiles,1),[]);
n_sm_pars = size(profiles,1);
threshold = chi2inv(0.95,n_sm_pars);
indices = zeros(n_sm_pars,size(profiles,2));
for i = 1:n_sm_pars
    for j = 1:size(profiles,2)
        profile = profiles{i,j};
        min_val = min(profile(end,:));
        indices(i,j) = sum(profile(end,[1,end]) >= min_val + threshold);
        % figureOnRight
        % for k = 1:4
        %     nexttile
        %     plot(profile(i,:), profile(k,:))
        % end
        % yline(min_val + threshold)
        % title(indices(i,j))
    end
end

