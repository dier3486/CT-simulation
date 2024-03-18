function [nextdata, writenum] = pooldatacopy(currdata, nextdata, ReadPoint, WritePoint, writenum, poolfields, force_flag)
% to copy the data from currunt data to next data in pipeline.
%   [nextdata, writenum] = pooldatacopy(currdata, nextdata, ReadPoint, WritePoint, writenum, poolfields, force_flag);
% or
%   [nextdata, writenum] = pooldatacopy(currdata, nextdata, ReadPoint, WritePoint, writenum);
% where the ReadPoint = status.pipepool.(currnode).ReadPoint, WritePoint = status.pipepool.(nextpool).WritePoint, to copy the
% data from input pool to output pool. Or to set the ReadPoint to your private 'currpool' reading position to copy the data
% from the private 'currpool' to output pool.
% Only one set of arrays can be transfered between the nodes, that means when you transfer the 'rawdata' and 'rawhead' the 
% 'image' can not be transfered because they are in different number. Otherwise, when you transfer the 'image' the 'rawdata' 
% and 'rawhead' have to be abandoned.

if writenum<=0
    writenum = 0;
%     return % do not return! empty fields shall be written to nextpool
end
if nargin < 6  || isempty(poolfields)
    poolfields = fieldnames(currdata);
end
if nargin < 7
    force_flag = false;
end

for ii = 1:length(poolfields)
    if isempty(poolfields{ii})
        % pass the {''}
        continue
    end
    if ~isfield(currdata, poolfields{ii})
        % pass the not exist fields
        continue;
    end
    if isfield(nextdata, poolfields{ii})
        if size(currdata.(poolfields{ii}), 2) == 1 && isstruct(currdata.(poolfields{ii}))
            % to recurse
            nextdata.(poolfields{ii}) = pooldatacopy(currdata.(poolfields{ii}), nextdata.(poolfields{ii}), ...
                ReadPoint, WritePoint, writenum, {}, true);
            continue
        end
        if isempty(currdata.(poolfields{ii}))
            % pass the empty data
            continue
            % some times we need to bypass empty fields in pool
        end
        % copy data
        nextdata.(poolfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
            currdata.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1);
    elseif force_flag
        % force to write
        if size(currdata.(poolfields{ii}), 2) == 1 && isstruct(currdata.(poolfields{ii}))
            % to recurse
            nextdata.(poolfields{ii}) = struct();
            nextdata.(poolfields{ii}) = pooldatacopy(currdata.(poolfields{ii}), nextdata.(poolfields{ii}), ...
                ReadPoint, WritePoint, writenum, {}, true);
            continue;
        end
        1;
        nextdata.(poolfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
            currdata.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1);
    end
end

end