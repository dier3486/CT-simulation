function u = TVpenalty_test1(f0, mu, lambda, u0, Crange)

% Crange = 1000+[-100 100];
if nargin<5
    Crange = [-inf inf];
end
Niter = 100;
tol_iter = 1e-2;
imgsize = size(f0);

f1 = f0.*1000;

s1 = (f1>=Crange(1)) & (f1<=Crange(2));
N1 = sum(s1(:));
f1(f1<Crange(1)) = Crange(1);
f1(f1>Crange(2)) = Crange(2);

b = zeros([imgsize 2]);
d = zeros([imgsize 2]);
delta = zeros(1, Niter);

if nargin<4 || isempty(u0)
    u0 = f1;
else
    u0 = u0.*1000;
    u0(~s1) = f1(~s1);
    [d, b] = fundbyub(u0, b, lambda);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, b, d, mu, lambda);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol_iter
        break;
    end
    [d, b] = fundbyub(u, b, lambda);
end

if delta(end)>delta(end-1)
    error('TV failed!');
end

u = u./1000;
u(~s1) = f0(~s1);
% delta(ii)
disp([ii delta(ii)]);
end

function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)
% d_y^{k+1} = shrink(D_yu^{k+1}+b_y^{k}, 1/\lambda)

[nx, ny] = size(u);

du = zeros(nx, ny, 2);
du(:,:,1) = [u(2:end,:) - u(1:end-1, :); zeros(1, ny)];
du(:,:,2) = [u(:, 2:end) - u(:, 1:end-1), zeros(nx, 1)];

b0 = b0 + du;
d1 = max(abs(b0)-1/lambda, 0).*(b0./abs(b0));
d1 = fillmissing(d1, 'constant', 0);

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, lambda)
[nx, ny] = size(u);

u1 = u([2:nx nx], :) + u([1 1:nx-1], :) + u(:, [2:ny ny]) + u(:, [1 1:ny-1]);
u1 = u1 + d([1 1:nx-1], :, 1) - d(:, :, 1) + d(:, [1 1:ny-1], 2) - d(:,:,2);
u1 = u1 - b([1 1:nx-1], :, 1) + b(:, :, 1) - b(:, [1 1:ny-1], 2) + b(:,:,2);
u1 = u1.*(lambda./(mu+4*lambda)) + f.*(mu./(mu+4*lambda));

end