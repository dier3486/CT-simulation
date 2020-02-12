% x0 = phantom.*1000;
% theta = 0:0.5:179.5;
% b = radon(x0, theta);
% 
% Nv = size(theta(:),1);

m = 10; n = 8;
A = rand(m, n);
b = rand(m,1);
x0 = zeros(n,1);

N = 20;
x = zeros(n, N+1);
x(:,1) = x0;
r = zeros(n, N);
% err = nan(N,1);
a = 0.5;
H = eye(n);
for ii = 1:N
    r(:,ii) = A'*A*x(:,ii) - A'*b;
    if ii>1
        y = r(:,ii) - r(:,ii-1);
        H = H + (dx - H*y)*dx'*H./(dx'*H*y);
    end
    dx = -(H * r(:, ii)).*a;
    x(:,ii+1) = x(:,ii) + dx;
end

err = sqrt(sum(r.^2,1));
figure;
plot(err, '.-');