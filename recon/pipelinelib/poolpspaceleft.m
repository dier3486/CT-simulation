function d = poolpspaceleft(currpool, viewcommon, viewexpand)
% to check how many space left in a pool

if nargin<2 || isempty(viewcommon)
    viewcommon = 1;
end

if nargin<3 || isempty(viewexpand)
    viewexpand = 0;
end

if isempty(currpool)
    % a Null pool
    d = inf;
    % A null pool have infinite space.
elseif isfield(currpool, 'WriteStuck') && currpool.WriteStuck
    % writing locked
    d = 0;
elseif currpool.circulatemode
    % circulate
%     if currpool.WriteEnd<1 || currpool.WriteEnd>currpool.poolsize
%         % no WriteEnd
%         d = inf;
%         % Note: We claim the size of a circulated buffer free from the limitation of WriteEnd is infinite. So when we want the
%         % pool 'automatically' locked after a rotation of data is filled we must set suitable poolsize and WriteEnd.
%     else
%         % to WriteEnd
%         d = mod(currpool.WriteEnd - currpool.WritePoint, currpool.poolsize) + 1;
%     end
    d = currpool.WriteEnd - currpool.WritePoint + 1;
    d = max(d, 0);
else
    % non-circulate
    % to WriteEnd or poolsize
    d = min(currpool.WriteEnd, currpool.poolsize - viewexpand) - currpool.WritePoint + 1;
    % no negtive
    d = max(d, 0);
end

if nargin >= 2
    d = floor(d/viewcommon) * viewcommon;
end

end