function out = flattenStruct(S)

out = struct();
fn = fieldnames(S);
for i = 1:length(fn)
    if ~isstruct(S.(fn{i}))
        out.(fn{i}) = S.(fn{i});
    else
        out = diveIn(S.(fn{i}),fn{i},out);
    end
end

end

function out = diveIn(S,name,out)

fn = fieldnames(S);
for i = 1:length(fn)
    
    if ~isstruct(S.(fn{i}))
        out.([name,'_DOT_',fn{i}]) = S.(fn{i});
    else
        out = diveIn(S.(fn{i}),fn{i},out);
    end

end

end