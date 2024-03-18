function d = poolp2p(p1, p2, psize, circ_flag)
% to calculte d = p2 - p1 in a pool
% d = poolp2p(p1, p2, currpool.poolsize, currpool.circulatemode);

if circ_flag
    if p2<1 || p2>psize
        d = inf;
    else
        d = mod(p2-p1, psize);
    end
else
    d = p2 - p1;
end

end