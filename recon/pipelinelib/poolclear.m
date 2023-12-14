function currpool = poolclear(currpool, poolfields)
% to clear the pool
%   currpool = poolclear(currpool);
% or
%   currpool = poolclear(currpool, poolfields);

if nargin < 2  || isempty(poolfields)
    poolfields = fieldnames(currpool);
end

for ii = 1:length(poolfields)
    if isfield(currpool, poolfields{ii})
        if isstruct(currpool.(poolfields{ii}))
            if size(currpool.(poolfields{ii}), 2) == 1
                % to recurse
                currpool.(poolfields{ii}) = poolclear(currpool.(poolfields{ii}));
            else
                % empty structure
                currpool.(poolfields{ii}) = struct([]);
            end
        else
            % empty
            currpool.(poolfields{ii}) = cast([], 'like', currpool.(poolfields{ii}));
        end
    end
end

end