function mirrorpool = poolmirror(currpool, mirrorpool)
% copy the fields of a pool expect the data

if nargin < 2
    mirrorpool = struct();
end

poolfields = fieldnames(currpool);

for ii = 1 : length(poolfields)
    switch poolfields{ii}
        case {'datafields', 'data', 'trace'}
            % skip
        otherwise
            mirrorpool.(poolfields{ii}) = currpool.(poolfields{ii});
    end

end

end