function sm_data = solveSM(sm,p,tt,C,D,condition_on_previous,resample_t)

persistent sm_pre_processor;
persistent sm_solver;
persistent sm_post_processor;

if isempty(sm_pre_processor)
    if isfield(sm,"pre_processor")
        sm_pre_processor = sm.pre_processor;
    else
        sm_pre_processor = @identityFunction;
    end
end

if isempty(sm_solver)
    switch sm.type
        case "ode"
            if condition_on_previous
                sm_solver = @solveODEWithExtend;
            else
                sm_solver = @solveODE;
            end
            if ~isfield(sm,"post_processor")
                sm_post_processor = @(x) x;
            else
                sm_post_processor = sm.post_processor;
            end
    end

end

[sm,p,tt,C,~] = sm_pre_processor(sm,p,tt,C,D);
sm_data = sm_solver(sm,p,tt,C,resample_t);
sm_data = sm_post_processor(sm_data);
end

function sm_data = solveODE(sm,p,tt,C,resample_t)
if isempty(resample_t)
    sm_data = solveODEHelper(sm,p,tt,C);
else
    sm_data = solveODEHelper(sm,p,resample_t,C);
end
end

function sm_data = solveODEHelper(sm,p,tt,C)
persistent tvals;
if isempty(tvals)
    tvals = unique([0; tt(:)]);
end
[~,sm_data] = sm.solver(@(t,x) sm.fn(t,x,p,C),tvals,sm.y0);
switch length(tt)
    case 1
        sm_data = sm_data(end,:);
    case 2
        if tt(1)==0
            sm_data = sm_data([1,end],:);
        else
            sm_data = sm_data(2:3,:);
        end
    otherwise
        if tt(1)~=0
            sm_data = sm_data(2:end,:);
        end
end
end

function sm_data = solveODEWithExtend(sm,p,tt,C,D,resample_t)
% final argument is resample_t and can be included to make this faster
sm_data = sm.solver(@(t,x) sm.fn(t,x,p,C),[0,tt(1)],sm.y0);
for i = 2:length(tt)
    sm_data = odextend(sm_data, [], tt(i), D.A(i-1,:));
end
if isempty(resample_t)
    sm_data = deval(tt,sm_data)';
else
    sm_data = deval(resample_t,sm_data)';
end
end

function varargout = identityFunction(varargin)
    varargout = varargin;
end