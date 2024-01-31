function currpool = poolrecycle(currpool, removenumber, poolfields)
% to recycle a data pool to nextpool in pipeline.
%   currpool = poolrecycle(currpool, removenumber, poolfields);
% or
%   currpool = poolrecycle(currpool, removenumber);

if nargin < 3  || isempty(poolfields)
    poolfields = fieldnames(currpool);
end

if removenumber<0
    removenumber = 0;
end

for ii = 1:length(poolfields)
    if isempty(poolfields{ii})
        % pass the {''}
        continue
    end
    if isfield(currpool, poolfields{ii})
        if size(currpool.(poolfields{ii}), 2) == 1 && isstruct(currpool.(poolfields{ii}))
            % to recurse
            currpool.(poolfields{ii}) = poolrecycle(currpool.(poolfields{ii}), removenumber);
            continue
        end
        if isempty(currpool.(poolfields{ii}))
            % pass the empty data
            continue
            % some times we need to bypass empty fields in pool
        end
        % recycle
        movesize = size(currpool.(poolfields{ii}), 2) - removenumber;
        currpool.(poolfields{ii})(:, 1:movesize) = currpool.(poolfields{ii})(:, removenumber+1 : end);
        currpool.(poolfields{ii})(:, movesize+1:end) = 0;
    end
end

end