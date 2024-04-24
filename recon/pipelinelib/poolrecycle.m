function [currpool, currdata] = poolrecycle(currpool, currdata, flag_unclock)
% to recycle a data pool to nextpool in pipeline.

if nargin<2 || isempty(currdata)
    currdata = currpool.data;
    % if you want to skip clearing currpool.data you shall set the currdata=struct() but not [].
end
if nargin<3
    flag_unclock = true;
    % it will unlock and recyle a writing stucked pool when it was used out, be careful.
end

if isfield(currpool, 'datafields')
    datafields = currpool.datafields;
else
    datafields = {};
end

% recyle depending on recylestrategy
switch currpool.recylestrategy
    case 0
        % never recycle
        % do nothing
    case 1
        % always recycle
        % recycle current pool
        [currpool, currdata] = poolrecycle2(currpool, currdata, currpool.ReadPoint-1, datafields);
    case 2
        if currpool.WritePoint >= currpool.warningstage
            [currpool, currdata] = poolrecycle2(currpool, currdata, currpool.ReadPoint-1, datafields);
        end
    otherwise
        error('Illeagal recylestrategy %d!', currpool.recylestrategy);
end

% unlock and force to recyle a stucked pool
WriteStuck = isfield(currpool, 'WriteStuck') && currpool.WriteStuck;
if flag_unclock && WriteStuck && currpool.AvailNumber==0
    currpool.WriteStuck = false;
    % all data = 0;
    [~, currdata] = poolrecycle2([], currdata, currpool.poolsize);
end

end