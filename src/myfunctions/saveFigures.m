function saveFigures(f,opts)

arguments
    f
    opts.save_figs = true; % if this function is called, then save figures by default
    opts.reprint = false; % whether or not to reprint figures
    opts.reprint_warning = true; % whether or not to warn about not reprinting
    opts.file_types = "fig"; % default to only saving a MATLAB figure file
    opts.nameless_warning = true; % whether to print warning for a nameless figure that will not be saved
    opts.fig_names = strings(0,1); % allow fig names to be input as options string
    opts.subfolder = "."; % subfolder within, e.g., figures/fig to hold files

    opts.resolution = '-r0'; % use screen resolution as default



end
% opts = defaultSaveFiguresOptions;
% if nargin>1 && ~isempty(input_opts)
%     opts = overrideDefaultOptions(opts,input_opts);
% end

if ~opts.save_figs
    return;
end

n_file_types = length(opts.file_types);
fig_folders = strings(n_file_types,1);
for j = 1:n_file_types
    fig_folders(j) = sprintf("figures/%s/%s",opts.file_types(j),opts.subfolder);
    if ~exist(fig_folders(j),"dir")
        mkdir(fig_folders(j))
    end
end

for i = 1:numel(f)
    if ~ishandle(f(i))
        continue; % do not save an empty graphics handle
    end

    if isempty(f(i).Name)
        % Only print named figs, and this fig is not named
        if numel(opts.fig_names)<i || isempty(opts.fig_names(i)) % nor was a name supplied in input_opts
            if opts.nameless_warning
                warning("Did not print Figure %d because no name was supplied.\n",f(i).Number)
            end
            continue
        else
            f(i).Name = opts.fig_names(i);
        end
    end

    for j = 1:n_file_types
        file_name = sprintf("%s/%s.%s",fig_folders(j),f(i).Name,opts.file_types(j));
        switch opts.file_types(j)
            case "fig"
                if opts.reprint || ~exist(file_name,"file")
                    savefig(f(i),file_name)
                elseif opts.reprint_warning
                    warning("Did not reprint %s because it already existed.\n" + ...
                        "  Set input_opts.reprint=true to force a reprint.\n" + ...
                        "  Turn this warning off by setting input_opts.reprint_warning=false.\n",file_name)
                end
            otherwise
                if opts.reprint || ~exist(file_name,"file")
                    print(opts.resolution,f(i),file_name,sprintf("-d%s",opts.file_types(j)))
                elseif opts.reprint_warning
                    warning("Did not reprint %s because it already existed.\n" + ...
                        "  Set input_opts.reprint=true to force a reprint.\n" + ...
                        "  Turn this warning off by setting input_opts.reprint_warning=false.\n",file_name) %#ok<*SPWRN> % this suppresses the warning about `warning` accepting `sprintf` inputs
                end
        end
    end

end
end

% function default_options = defaultSaveFiguresOptions
% 
% default_options.save_figs = true; % if this function is called, then save figures by default
% default_options.reprint = false; % whether or not to reprint figures
% default_options.reprint_warning = true; % whether or not to warn about not reprinting
% default_options.file_types = "fig"; % default to only saving a MATLAB figure file
% default_options.nameless_warning = true; % whether to print warning for a nameless figure that will not be saved
% default_options.fig_names = strings(0,1); % allow fig names to be input as options string
% default_options.subfolder = "."; % subfolder within, e.g., figures/fig to hold files
% 
% default_options.resolution = '-r0'; % use screen resolution as default
% 
% end

