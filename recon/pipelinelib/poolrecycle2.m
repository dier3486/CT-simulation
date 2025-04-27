function [currpool, currdata, removenumber] = poolrecycle2(currpool, currdata, removenumber, poolfields)
% to recycle a data pool to nextpool in pipeline.
%   [currpool, currdata] = poolrecycle2(currpool, currdata, removenumber, poolfields);

% default removenumber
if nargin<3 || isempty(removenumber)
    removenumber = currpool.ReadPoint - currpool.keepbottom - 1;
end
if removenumber<0
    removenumber = 0;
end

% default fields
if nargin < 4 || isempty(poolfields)
    poolfields = fieldnames(currdata);
end

for ii = 1:length(poolfields)
    if isempty(poolfields{ii})
        % pass the {''}
        continue
    end
    if isfield(currdata, poolfields{ii})
        if size(currdata.(poolfields{ii}), 2) == 1 && isstruct(currdata.(poolfields{ii}))
            % to recurse
            [~, currdata.(poolfields{ii})] = poolrecycle2([], currdata.(poolfields{ii}), removenumber);
            continue
        end
        if isempty(currdata.(poolfields{ii}))
            % pass the empty data
            continue
            % some times we need to bypass empty fields in pool
        end
        % recycle
        movesize = size(currdata.(poolfields{ii}), 2) - removenumber;
        if movesize > 0
            currdata.(poolfields{ii})(:, 1:movesize) = currdata.(poolfields{ii})(:, removenumber+1 : end);
            currdata.(poolfields{ii})(:, movesize+1:end) = 0;
        else
            currdata.(poolfields{ii})(:, :) = 0;
        end
    end
end

% reset the currpool
if ~isempty(currpool)
    % move ReadPoints
    currpool.ReadPoint = currpool.ReadPoint - removenumber;
    currpool.ReadStart = currpool.ReadStart - removenumber;
    currpool.ReadEnd = currpool.ReadEnd - removenumber;
    % move WritePoints
    currpool.WritePoint = currpool.WritePoint - removenumber;
    currpool.WriteStart = currpool.WriteStart - removenumber;
    currpool.WriteEnd = currpool.WriteEnd - removenumber;
    % The ReadViewindex and AvailViewindex is not changed.
    % move AvailPoint
    currpool.AvailPoint = currpool.AvailPoint - removenumber;
end
% The currpool can be empty. If you want to skip resetting the points in currpool just let the input in empty.

end