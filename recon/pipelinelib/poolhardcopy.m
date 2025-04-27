function outdata = poolhardcopy(outdata, indata, index_out, index_in, copyfields, method)

if nargin < 5 || isempty(copyfields)
    copyfields = fieldnames(indata);
end
if nargin < 6 || isempty(method)
    method = 'overwrite';
    % or 'cum'
end

for ii = 1:length(copyfields)
    if isempty(copyfields{ii})
        % pass the {''}
        continue
    end
    if ~isfield(indata, copyfields{ii})
        % pass the not exist fields
        continue;
    end
    if size(indata.(copyfields{ii}), 2) == 1 && isstruct(indata.(copyfields{ii}))
        % recurse the struct
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = struct();
        end
        outdata.(copyfields{ii}) = poolhardcopy(outdata.(copyfields{ii}), indata.(copyfields{ii}), index_out, index_in, [], method);
    else
        if ~isfield(outdata, copyfields{ii}) || isempty(outdata.(copyfields{ii}))
            outdata.(copyfields{ii}) = cast([], 'like', indata.(copyfields{ii}));
            outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
        else
            if size(outdata.(copyfields{ii}), 2) < max(index_out)
                outdata.(copyfields{ii})(:, max(index_out)) = 0;
            end
            switch lower(method)
                case 'overwrite'
                    outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
                case 'cum'
                    outdata.(copyfields{ii})(:, index_out) = ...
                        outdata.(copyfields{ii})(:, index_out) + indata.(copyfields{ii})(:, index_in);
                case 'clear'
                    outdata.(copyfields{ii})(:, index_out) = 0;
                otherwise
                    % do nothing
                    0;
                    % We can support function handles here, later.
            end
        end
        % note: the struct vector is not fully supported.
    end
end

end