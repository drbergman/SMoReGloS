function profiles = cleanProfiles(files,options)

% this function will take in a collection of profiles (out) along with a
% threshold and "clean" them. Clean means make sure the value does not
% cross the threshold more than once in either direction AND trim the ends
% so that we only have the 95% CI and not however far over that we stepped.

arguments
    files struct
    options.boundary_tolerance double = 0
end

load(files.profiles,"profiles","profile_params") % profiles from ABM
threshold = chi2inv(0.95,size(profiles,1));

sz = size(profiles);
npars = size(profiles,1);
profiles = reshape(profiles,npars,[]);
n_abm_vecs = size(profiles,2);

for i = 1:n_abm_vecs
    for j = 1:npars
        profiles{j,i} = cleanThisProfile(profiles{j,i},i,j,threshold);
        if any(options.boundary_tolerance > 0)
            profiles{j,i} = applyBoundaryThreshold(j,profiles{j,i},profile_params,options.boundary_tolerance);
        end
    end
end

profiles = reshape(profiles,sz);

end

function profile = cleanThisProfile(profile,i,j,threshold)

try % if this returns without error, then just trim the profile
    profile = trimProfile(profile,threshold);
catch % but if there's an error, then there's multiple crosses, remove extraneous times it rose above the threshold
    [inds_before_drop,inds_before_rise,max_val] = findMaxCrosses(profile(end,:),threshold);
    [i1,~] = find(inds_before_drop'-inds_before_rise==1); % I1 contains indices where the value began to drop; if these come right after an increase over the threshold, then these are likely outliers we want to ignore
    if ~isempty(i1) % then there was at least one time a jump occurred, ignore these
        profile(:,inds_before_drop(i1)) = [];
        profile = cleanThisProfile(profile,i,j,threshold);
    else
        % plot unclean profile
        f=figureOnRight;
        plot(profile(j,:),profile(end,:),"Marker",".")
        yline(max_val)

        [ST,current_ws_ind] = dbstack();
        stopText = strsplit(sprintf('in %s at %d', ST(current_ws_ind).name, ST(current_ws_ind).line + 3));
        dbstop(stopText{:});
        fprintf("Unknown error encountered in trimming profile %d at parameter %d. \n   You can fix manually now and then continue.\n",i,j)

        if false %#ok<*UNRCH> % suppress warning that the following lines cannot be reached, they are here to make them easily copy-pasted into Command Window
            % if just multiple one-step peaks over max_val (this can no longer be the fix since I have changed the check above from numel(i1)==1 to ~isempty(i1), meaning I will remove all one-step peaks)
            profile(:,inds_before_drop(i1)) = [];

            % if there's a single, multi-step peak over max_val
            profile(:,(inds_before_rise(i2_ind)+1):inds_before_drop(i1_ind)) = [];
            
        end
        hold on;
        plot(profile(j,:),profile(end,:))
        profile = cleanThisProfile(profile,i,j,threshold);
        close(f)
    end
end

end

function profile = applyBoundaryThreshold(profile_par_index,profile,profile_params,boundary_tolerance)

lb = profile_params.lb;
ub = profile_params.ub;
if size(boundary_tolerance,2)==2
    tol = boundary_tolerance;
else
    d = ub-lb;
    tol = boundary_tolerance .* d;
    tol = [tol,tol];
end

for i = 1:(size(profile,1)-1)
    if i==profile_par_index
        continue
    end
    % idx = profile(i,:) <= lb(i) + tol(i,1) | profile(i,:) >= ub(i) - tol(i,2);
    % if any(idx)
    %     figureOnRight();
    %     for k = 1:size(profile,1)
    %         subplot(size(profile,1),1,k)
    %         hold on;
    %         plot(profile(profile_par_index,~idx),profile(k,~idx),"Marker","+")
    %         plot(profile(profile_par_index,idx),profile(k,idx),"Marker","+")
    %     end
    %     disp('')
    % end

    % check if any parameter approaches boundary within tolerance at end points
    while true
        if profile(end,1) > min(profile(end,:)) && (profile(i,1) <= lb(i) + tol(i,1) || profile(i,1) >= ub(i) - tol(i,2))
            profile = profile(:,2:end);
        else
            break
        end
    end
    while true
        if profile(end,end) > min(profile(end,:)) && (profile(i,end) <= lb(i) + tol(i,1) || profile(i,end) >= ub(i) - tol(i,2))
            profile = profile(:,1:end-1);
        else
            break
        end
    end
end
end
