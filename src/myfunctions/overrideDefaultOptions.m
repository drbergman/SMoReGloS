function out = overrideDefaultOptions(default_options,new_options)

% overrides the default options structure default_options with those in
% options. If any fields of options are not included in default_options,
% then this returns an error.

if ~isequal(class(default_options),class(new_options))
    error("Options are not in the same class.")
end

switch class(default_options)
    case 'struct'
        fn_default = fieldnames(default_options);
        fn_new = fieldnames(new_options);

        unexpected_fields = setdiff(fn_new,fn_default);
        if ~isempty(unexpected_fields)
            fprintf("The following fields are included that are not found in the default options:\n")
            for i = 1:numel(unexpected_fields)
                fprintf("  %d. %s\n",i,unexpected_fields{i})
            end
            error("Unknown fields in input.")
        end

        out = default_options;

        for i = 1:numel(fn_new)
            out.(fn_new{i}) = new_options.(fn_new{i});
        end

    case 'containers.Map'
        keys_default = default_options.keys;
        keys_new = new_options.keys;

        unexpected_keys = setdiff(keys_new,keys_default);
        if ~isempty(unexpected_keys)
            fprintf("The following keys are included that are not found in the default options:\n")
            for i = 1:numel(unexpected_keys)
                fprintf("  %d. %s\n",i,unexpected_keys{i})
            end
            error("Unknown keys in input.")
        end

        out = containers.Map(default_options.keys,default_options.values);

        for i = 1:numel(keys_new)
            out(keys_new{i}) = new_options(keys_new{i});
        end
        
    otherwise
        error("Unsure how to handle options as class %s.",class(default_options))

end