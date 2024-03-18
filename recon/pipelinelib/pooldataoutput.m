function [currpool, nextpool, Rinfo] = pooldataoutput(currpool, nextpool, viewcommon)
% to copy the data from currpool to nextpool
% a simplified version in using pooldatacopy
% [currpool, currdata, nextpool, nextdata, output] = pooldataoutput(currpool, currdata, nextpool, nextdata);

if nargin < 2
    viewcommon = 1;
end

% call pooldataoutput2
[currpool, nextpool, nextpool.data, Rinfo] = pooldataoutput2(currpool, currpool.data, nextpool, nextpool.data, ...
    viewcommon);

end