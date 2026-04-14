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
            % shall not happen on deployment
            outdata.(copyfields{ii}) = zeros(size(indata.(copyfields{ii}), 1), 0, 'like', indata.(copyfields{ii}));
            outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
            % pass the lib.pointers
        elseif isa(outdata.(copyfields{ii}), 'lib.pointer') || isa(indata.(copyfields{ii}), 'lib.pointer')
            continue;
        else
            if size(outdata.(copyfields{ii}), 2) < max(index_out)
                outdata.(copyfields{ii})(:, max(index_out)) = 0;
            end
            % while the outdata isreal only copy the real part
            outreal = isreal(outdata.(copyfields{ii}));
            inreal = isreal(indata.(copyfields{ii}));
            switch lower(method)
                case 'overwrite'
                    if outreal && ~inreal
                        outdata.(copyfields{ii})(:, index_out) = real(indata.(copyfields{ii})(:, index_in));
                    else
                        outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
                    end
                case 'cum'
                    if outreal && ~inreal
                        outdata.(copyfields{ii})(:, index_out) = ...
                            outdata.(copyfields{ii})(:, index_out) + real(indata.(copyfields{ii})(:, index_in));
                    else
                        outdata.(copyfields{ii})(:, index_out) = ...
                            outdata.(copyfields{ii})(:, index_out) + indata.(copyfields{ii})(:, index_in);
                    end
                case 'clear'
                    outdata.(copyfields{ii})(:, index_out) = 0;
                otherwise
                    % do nothing
                    0;
                    % We can support function handles here, later.
            end
            % to fix the auto complex -> real bug (need not in deployment codes)
            if ~outreal && isreal(outdata.(copyfields{ii}))
                % it is a bug of matlab, I can not even avoid it by using cast or more if-else.
                outdata.(copyfields{ii}) = complex(outdata.(copyfields{ii}));
            end
            % I know the GPU/CPU can be autoly gathered or put.
        end
        % note: the struct vector is not fully supported.
    end
end

end