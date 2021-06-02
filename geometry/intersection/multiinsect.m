function D = multiinsect(L, R)

D = min(R, [], 2) - max(L, [], 2);
D(D<0) = 0;
