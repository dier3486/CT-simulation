function currdata = poolresize(currdata, sizeto, resizefields, refdata, force_flag)
% resize the buffer

if nargin < 4 || isempty(refdata)
    refdata = currdata;
end

if nargin < 3 || isempty(resizefields)
    resizefields = fieldnames(refdata);
end

if nargin < 5
    force_flag = true;
end

for ii = 1:length(resizefields)
    if isempty(resizefields{ii})
        % pass the {''}
        continue
    end
    if ~isfield(currdata, resizefields{ii}) && ~isfield(refdata, resizefields{ii})
        % pass the not exist fields
        continue;
    end
    if ~isfield(currdata, resizefields{ii}) || (isempty(currdata.(resizefields{ii})) && ...
            size(currdata.(resizefields{ii}), 1) == 0)
        if isfield(refdata, resizefields{ii})
            if ~isstruct(refdata.(resizefields{ii}))
                currdata.(resizefields{ii}) = refdata.(resizefields{ii})(:, []);
            else
                currdata.(resizefields{ii}) = struct();
            end
        else
            % ??
        end
    end
    if size(currdata.(resizefields{ii}), 2) == 1 && isstruct(currdata.(resizefields{ii}))
        % to recurse
        if isfield(refdata, resizefields{ii})
            currdata.(resizefields{ii}) = poolresize(currdata.(resizefields{ii}), sizeto, [], ...
                refdata.(resizefields{ii}), force_flag);
        else
            currdata.(resizefields{ii}) = poolresize(currdata.(resizefields{ii}), sizeto, [], [], force_flag);
        end
        continue;
    end
    % resize
    if size(currdata.(resizefields{ii}), 2) < sizeto
        currdata.(resizefields{ii})(:, sizeto) = 0;
    elseif size(currdata.(resizefields{ii}), 2) > sizeto && force_flag
        currdata.(resizefields{ii})(:, sizeto+1:end) = [];
    end
    

end

end