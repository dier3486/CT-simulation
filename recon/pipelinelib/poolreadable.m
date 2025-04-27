function d = poolreadable(currpool, minlimit)
% to check how many data readable in a pool

if nargin<2
    % readable number shall >= minlimit, or retrun 0
    minlimit = 1;
end

if isempty(currpool)
    % a Null pool
    d = 0;
elseif isfield(currpool, 'ReadStuck') && currpool.ReadStuck
    % reading locked
    d = 0;
else
%     RP2RE = poolp2p(currpool.ReadPoint, currpool.ReadEnd, currpool.poolsize, currpool.circulatemode);
%     d = min(currpool.AvailNumber, RP2RE+1);
%     d = floor(d/viewcommon) * viewcommon;
    d = min(currpool.ReadEnd, currpool.AvailPoint) - currpool.ReadPoint + 1;
    if d < min(minlimit, 0)
        d = 0;
    end
    % In a buffer the readable data length (d) is the AvailNumber starting from the ReadPoint while AvailNumber>=minlimit, or
    % it is 0.
    % But it is still possible for a node to visit the data out of that range especialy while a buffer is in circulate mode.
end

% Note: This function did not consider to fit the d into viewcommon, user can do that by this way,
%   d = floor(d/viewcommon) * viewcommon;
% There is very strong reason for why we don't employ the viewcommon in this function, do not easily change that.

end