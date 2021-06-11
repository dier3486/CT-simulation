function L = splaplace2D(m, n)
% 2D sparse Laplace, m*n

if nargin<2
    n = m;
end

% Lx
Lx1 = cell(1, n);
Lx1(:) = {spdiags(ones(m-1,1)*[-1 1], [0 1], m-1, m)};
Lx = blkdiag(Lx1{:});

% Ly
Ly1 = cell(1, m);
Ly1(:) = {spdiags(ones(n-1,1)*[-1 1], [0 1], n-1, n)};
Ly = blkdiag(Ly1{:});
index_y(reshape(1:m*n, m, n)') = 1:m*n;
Ly = Ly(:, index_y(:));

L = Lx'*Lx + Ly'*Ly;
end