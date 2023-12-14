function [nextpool, writenum] = pooldatacopy(currpool, nextpool, ReadPoint, WritePoint, writenum, poolfields, force_flag)
% to copy the data from currpool to nextpool in pipeline.
%   [nextpool, writenum] = pooldatacopy(currpool, nextpool, ReadPoint, WritePoint, writenum, poolfields, force_flag);
% or
%   [nextpool, writenum] = pooldatacopy(currpool, nextpool, ReadPoint, WritePoint, writenum);
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
    poolfields = fieldnames(currpool);
end
if nargin < 7
    force_flag = false;
end

for ii = 1:length(poolfields)
    if isempty(poolfields{ii})
        % pass the {''}
        continue
    end
    if isfield(nextpool, poolfields{ii})
        if size(currpool.(poolfields{ii}), 2) == 1 && isstruct(currpool.(poolfields{ii}))
            % to recurse
            nextpool.(poolfields{ii}) = pooldatacopy(currpool.(poolfields{ii}), nextpool.(poolfields{ii}), ...
                ReadPoint, WritePoint, writenum, {}, true);
            continue
        end
        if isempty(currpool.(poolfields{ii}))
            % pass the empty data
            continue
            % some times we need to bypass empty fields in pool
        end
        % copy data
        nextpool.(poolfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
            currpool.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1);
%         poollength = size(nextpool.(poolfields{ii}), 2);
%         if poollength >= (WritePoint + writenum - 1)
%             nextpool.(poolfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
%                 currpool.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1);
%         elsei
%             nextpool.(poolfields{ii}) = [nextpool.(poolfields{ii})(:, 1:WritePoint-1) ...
%                 currpool.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1)];
%         end
    elseif force_flag
        % force to write
        if size(currpool.(poolfields{ii}), 2) == 1 && isstruct(currpool.(poolfields{ii}))
            % to recurse
            nextpool.(poolfields{ii}) = struct();
            nextpool.(poolfields{ii}) = pooldatacopy(currpool.(poolfields{ii}), nextpool.(poolfields{ii}), ...
                ReadPoint, WritePoint, writenum, {}, true);
            continue;
        end
        1;
        nextpool.(poolfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
            currpool.(poolfields{ii})(:, ReadPoint : ReadPoint + writenum - 1);
    end
end

end