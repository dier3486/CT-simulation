function u = GlogTV1D(img0, mu, lambda, logC, u0, DIM, Crange, Niter, tol)
% TV with log G fun
% u = GlogTV1D(img0, mu, lambda, tanhC, u0, DIM, Crange, Niter, tol);
% or u = GlogTV1D(img0, mu, lambda, tanhC);

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

% f1 = img0;
if any(isinf(Crange))
    Crange(isinf(Crange)) = 1e8.*sign(Crange(isinf(Crange)));
    f1 = hardentanh(img0, 1e-6, Crange);
else
    C0 = sum(Crange)/2;
    f1 = hardentanh(img0 - C0, 1e-6, Crange(2)-C0) + C0;
end

N1 = prod(imgsize);

b = zeros(size(img0), 'like', img0);
d = zeros(size(img0), 'like', img0);
delta = zeros(1, Niter, 'like', img0);

if nargin<5 || isempty(u0)
    u0 = f1;
else
    u0 = permute(u0, dimindex);
    u0 = reshape(u0, imgsize(1), []);
    u0 = img0 + u0 - f1;
%     [d, b] = fundbyub(u0, b, lambda);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, b, d, mu, lambda, logC);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(u, b, lambda);
%     disp([ii delta(ii)]);
end

if delta(end)>delta(end-1)
    error('TV failed!');
end

u = img0 + u - f1;

% disp([ii delta(ii)]);

u = ipermute(reshape(u, imgsize), dimindex);

end


function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)

% boundary 1
du = b0.*0;
du(1:end-1, :) = u(2:end, :) - u(1:end-1, :);

% % boundary 2
% du = u([2:end end], :) - u([1:end-1, end-1], :);
% du(end, :) = -du(end, :);

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);
d1 = max(absb0-1/lambda, 0).*(b0./absb0);
% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function [u1, f1] = funG(f, u, b, d, mu, lambda, logC)

% boundary 1
u1 = u([2:end end], :) + u([1 1:end-1], :);
d = d - d([1 1:end-1], :);
b = b - b([1 1:end-1], :);

% % boundary 2
% u1 = u([2:end end-1], :) + u([2 1:end-1], :);
% 
% d(2:end,:) = d(2:end,:) - d(1:end-1,:);
% d(1,:) = d(1,:).*2;
% 
% b(2:end,:) = b(2:end,:) - b(1:end-1,:);
% b(1, :) = b(1, :).*2;

u1 = u1 - d + b;

% f1 = tanh((f-u)./tanhC).*tanhC + u;
% f1 = hardentanh2(f-u, 1.0, tanhC) + u;

f1 = log(abs(f-u)./logC+1).*logC.*sign(f-u) + u;

u1 = u1.*(lambda./(mu+2*lambda)) + f1.*(mu./(mu+2*lambda));

end