function out = r2c(nr,nc,i)

% take in number of rows and columns for subplots and the desired index
% using normal matlab indexing (first vary rows then columns) and output
% the subplot index for matlab's subplot indexing

switch numel(i)
    case 1 % input is given in normal matlab indexing
        if i>nr*nc
            error('index is too big')
        end
        [ri,ci] = ind2sub([nr,nc],i);
        out = sub2ind([nc,nr],ci,ri);
    case 2
        if i(1)>nr
            error('row index is too big')
        elseif i(2)>nc
            error('column index is too big')
        end
        out = i(2) + (i(1)-1)*nc;        
    otherwise
        error('need one or two scalars for this')
end