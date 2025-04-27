function [currpool, currdata, removenumber] = poolrecycle(currpool, currdata, removenumber, flag_unclock)
% to recycle a data pool to nextpool in pipeline.
%   [currpool, currdata, removenumber] = poolrecycle(currpool, currdata, removenumber, flag_unclock);
% or simplified, currpool = poolrecycle(currpool);

if nargin<4
    flag_unclock = true;
    % it will unlock and recyle a writing stucked pool when it was used out, be careful.
end

if length(currpool) > 1
    % simplified using of multi pools (do not employ it on carried pools)
    for ii = 1:length(currpool)
        currpool(ii) = poolrecycle(currpool(ii), [], [], flag_unclock);
    end
    return;
end

if nargin<2 || isempty(currdata)
    currdata = currpool.data;
    % if you want to skip clearing currpool.data you shall set the input currdata=struct() but not [], like this,
    % currpool = poolrecycle(currpool, struct(), removenumber, flag_unclock);
end

% view number to be removed
if nargin<3
    removenumber = [];
    % I know while it is empty the removenumber = currpool.ReadPoint-1.
end

if isfield(currpool, 'datafields')
    datafields = currpool.datafields;
else
    datafields = {};
end

% recyle depending on recylestrategy
if ~currpool.circulatemode
    switch currpool.recylestrategy
        case 0
            % never recycle
            removenumber = 0;
            % do nothing
        case 1
            % always recycle
            % recycle current pool
            [currpool, currdata, removenumber] = poolrecycle2(currpool, currdata, removenumber, datafields);
        case 2
            if currpool.WritePoint >= currpool.warningstage
                [currpool, currdata, removenumber] = poolrecycle2(currpool, currdata, removenumber, datafields);
            end
        otherwise
            error('Illeagal recylestrategy %d!', currpool.recylestrategy);
    end
else
    % for circulatemode, set the removenumber to 0 first.
    removenumber = 0;
end

if flag_unclock
    % unlock and force to recyle a stucked pool (mostly happen in the end of a shot)
    WriteStuck = isfield(currpool, 'WriteStuck') && currpool.WriteStuck;
    % writing locked and all available data be read
    AvailNumber = currpool.AvailPoint - currpool.ReadPoint + 1;
    if WriteStuck && AvailNumber==0
        % unlock writing
        currpool.WriteStuck = false;
        % remove read end
        if isavail(currpool.ReadEnd)
            currpool.ReadEnd = Inf;
        end
        % period draw back
        if currpool.circulatemode
            origReadPoint = currpool.ReadPoint;
            currpool.ReadPoint = mod(currpool.ReadPoint - 1, currpool.poolsize) + 1;
            currpool.WritePoint = currpool.WritePoint + currpool.ReadPoint - origReadPoint;
            % to return the removenumber
            removenumber = currpool.poolsize;
        end
        % reset read start
        currpool.ReadStart = currpool.ReadPoint;
        % remove write end (only for non-circulate mode)
        if isavail(currpool.WriteEnd) && ~currpool.circulatemode
            currpool.WriteEnd = Inf;
        end
        % reset write start
        currpool.WriteStart = currpool.WritePoint;
        % clear all data to 0
        [~, currdata] = poolrecycle2([], currdata, currpool.poolsize);
        % reset AvailPoint
        currpool.AvailPoint = currpool.ReadPoint - 1;
        % reset isshotstart
        currpool.isshotstart = true;
    end
    % Note: in principle a pool need to be re-initial before using after these 'unlock' operations.
end

if nargin<2 || isempty(currdata)
    % Note again, if you want to skip clearing currpool.data you shall set the currdata=struct() but not [].
    currpool.data = currdata;
end

end