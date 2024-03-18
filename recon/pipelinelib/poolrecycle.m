function [currpool, currdata] = poolrecycle(currpool, currdata)
% to recycle a data pool to nextpool in pipeline.

if nargin<2
    currdata = currpool.data;
end

if isfield(currpool, 'datafields')
    datafields = currpool.datafields;
else
    datafields = {};
end

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

% unlock
WriteStuck = isfield(currpool, 'WriteStuck') && currpool.WriteStuck;
if WriteStuck && currpool.AvailNumber==0
    currpool.WriteStuck = false;
end

end