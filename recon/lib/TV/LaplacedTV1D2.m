function u = LaplacedTV1D2(img0, mu, lambda, laplace, u0, DIM, Crange, Niter, tol)
% TV + Laplace, test2
% u = LaplacedTV1D(img0, mu, lambda, laplace, u0, DIM, Crange, Niter, tol)
% or u = LaplacedTV1D(f0, mu, lambda, laplace);
% typicaly, img0 is around 1000 HF, mu=0.1~1, lambda<mu/4, laplace=0~1

if nargin<6 || isempty(DIM)
    DIM = 1;
end
if nargin<7 || isempty(Crange)
    Crange = [-inf inf];
    % Crange = 1000+[-100 100];
end
if nargin<8 || isempty(Niter)
    Niter = 100;
end
if nargin<9 || isempty(tol)
    tol = 1e-2;
end

dimindex = 1 : ndims(img0);
dimindex = [DIM dimindex(dimindex~=DIM)];
img0 = permute(img0, dimindex);
% permuted size
imgsize = size(img0);
img0 = reshape(img0, imgsize(1), []);

f1 = img0;
s1 = (f1>=Crange(1)) & (f1<=Crange(2));
N1 = sum(s1(:));
f1(f1<Crange(1)) = Crange(1);
f1(f1>Crange(2)) = Crange(2);

b = zeros(size(img0), 'like', img0);
d = zeros(size(img0), 'like', img0);
delta = zeros(1, Niter, 'like', img0);

if nargin<5 || isempty(u0)
    u0 = f1;
else
    u0 = permute(u0, dimindex);
    u0 = reshape(u0, imgsize(1), []);
    u0(~s1) = f1(~s1);
    [d, b] = fundbyub(u0, b, lambda);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, b, d, mu, lambda, laplace);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(u, b, lambda);
end

if delta(end)>delta(end-1)
    error('TV failed!');
end

u(~s1) = img0(~s1);

% disp([ii delta(ii)]);

u = ipermute(reshape(u, imgsize), dimindex);

end


function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)

du = b0.*0;
du(1:end-1, :) = u(2:end, :) - u(1:end-1, :);

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);
d1 = max(absb0-1/lambda, 0).*(b0./absb0);
% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, lambda, laplace)
nx = size(u, 1);


u1 = u([2:nx nx], :) + u([1 1:nx-1], :);

L = laplace/mu;
u1 = u1.*(1 + L/2) - u.*L;

d1 = d - d([1 1:nx-1], :);
d1(1, :) = -d1(2, :)./2;
% d1(2:end, :) = d1(2:end, :) + d1(1:end-1, :);
% d1 = d1./2;

u1 = u1 - d1;

b1 = b - b([1 1:nx-1], :);
% b1(2:end, :) = b1(2:end, :) + b1(1:end-1, :);
% b1 = b1./2;
b1(1, :) = -b1(2, :)./2;

u1 = u1 + b1;


u1 = u1.*(lambda./(mu+2*lambda)) + f.*(mu./(mu+2*lambda));

end