function S = structcellcurse(S, curseindex)
% recurse the cell2mat() on the fields of structure S
% S = structcellcurse(S);
% or S = structcellcurse(S, n); n is the depth of the recurse

if nargin<2
    curseindex = inf;
end

if curseindex<0
    return
end

if iscell(S) && all(cellfun(@isstruct, S))
    try
        % try cell2mat()
        S = [S{:}];
    catch
        % do nothing
        1;
    end
end

if isstruct(S)
    if size(S(:), 1)>1
        % loop the structure list to recurse
        for ii = 1:size(S(:))
            S(ii) = structcellcurse(S(ii), curseindex);
        end
    else
        % loop the fields of S to recurse, curseindex-1
        for ifld = fieldnames(S)'
            S.(ifld{1}) = structcellcurse(S.(ifld{1}), curseindex-1);
        end
    end
end

if iscell(S)
    % loop the cell to recurse
    for ii = 1:size(S(:), 1)
        if isstruct(S{ii})
            % recurse a structure in cell
            S{ii} = structcellcurse(S{ii}, curseindex);
        elseif iscell(S{ii})
            % recurse a cell in cell
            S{ii} = structcellcurse(S{ii}, curseindex);
        end
    end
end

end
