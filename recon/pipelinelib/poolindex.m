function r = poolindex(currpool, index)
% to return the index in a pool

r = index(1) : index(2);

if currpool.circulatemode
    r = mod(r - 1, currpool.poolsize) + 1;
end


end