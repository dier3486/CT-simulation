function p = panalty1(A, alpha, beta)

[m,n] = size(A);

Dlr = A(:, 2:end) - A(:, 1:end-1);
Dud = A(2:end, :) - A(1:end-1, :);

p = zeros(m, n);
p(:, 2:end) = p(:, 2:end) + TVexp(Dlr, alpha);
p(:, 1:end-1) = p(:, 1:end-1) + TVexp(-Dlr, alpha);
p(2:end, :) = p(2:end, :) + TVexp(Dud, alpha);
p(1:end-1, :) = p(1:end-1, :) + TVexp(-Dud, alpha);

p = p.*beta;