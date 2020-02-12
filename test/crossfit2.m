function r = crossfit2(p, x, y, s, lambda)

[m, ~] = size(x);

d1 = 1-[0; p(:)]-[p(:); 0];
M = spdiags([[p; 0] d1 [0; p]], [-1 0 1], m, m);
r1 = x + log2(M\(2.^(-y)));

r = [r1(s); p(:).*lambda]';

end