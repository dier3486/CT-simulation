function [u, delta] = BregmanTV2D(img0, mu, Cl, u0, Crange, Niter, tol)
% Split Bregman method for TV
% u = BregmanTV2D(img0, mu, Cl, u0, Crange, Niter, tol)
% or u = BregmanTV2D(f0, mu, Cl);
% the Cl is lambda./mu
% typicaly, img0 is around 1000 HF, mu=0.1~1, lambda<mu/4


if nargin<5 || isempty(Crange)
    Crange = [-inf inf];
    % Crange = 1000+[-100 100];
end
if nargin<6 || isempty(Niter)
    Niter = 100;
end
if nargin<7 || isempty(tol)
    tol = 1e-2;
end
imgsize = size(img0);
if length(imgsize) < 3
    imgsize = [imgsize 1];
end 

f1 = img0;
s1 = (f1>=Crange(1)) & (f1<=Crange(2));
N1 = sum(s1(:));
f1(f1<Crange(1)) = Crange(1);
f1(f1>Crange(2)) = Crange(2);

GPUonoff = isa(img0, 'gpuArray');
fclass = classGPU(img0);

if GPUonoff
    b = zeros([imgsize 2], fclass, 'gpuArray');
    d = zeros([imgsize 2], fclass, 'gpuArray');
    delta = zeros(1, Niter, fclass, 'gpuArray');
else
    b = zeros([imgsize 2], fclass);
    d = zeros([imgsize 2], fclass);
    delta = zeros(1, Niter, fclass);
end

if nargin<4 || isempty(u0)
    u0 = f1;
else
    u0(~s1) = f1(~s1);
%     [d, b] = fundbyub(u0, b, lambda);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, b, d, Cl);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(u, b, mu.*Cl);
end

delta = delta(1:ii);

if delta(end)>delta(end-1)
    error('TV failed!');
end

% u = u./1000;
u(~s1) = img0(~s1);
% delta(ii)
% disp([ii delta(ii)]);
end


function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)
% d_y^{k+1} = shrink(D_yu^{k+1}+b_y^{k}, 1/\lambda)

[nx, ny, ~] = size(u);

% du = zeros(nx, ny, nz, 2);
du = b0.*0;
% du(1:end-1,:,:,1) = u(2:end,:,:) - u(1:end-1, :,:);
% du(:,1:end-1,:,2) = u(:, 2:end,:) - u(:, 1:end-1,:);

du(:,:,:, 1) = u([2:nx  nx], :,:) - u([1:nx-1  nx-1], :,:);
du(:,:,:, 2) = u(:, [2:ny  ny],:) - u(:, [1:ny-1  ny-1],:);

du(end,:,:, 1) = -du(end,:,:, 1);
du(:,end,:, 2) = -du(:,end,:, 2);

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);
d1 = max(absb0-1./lambda, 0).*(b0./absb0);
% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function u = funG(f, u, b, d, Cl)
[nx, ny, ~] = size(u);

u = u([2:nx nx], :, :) + u([1 1:nx-1], :, :) + u(:, [2:ny ny], :) + u(:, [1 1:ny-1], :);

% u1 = u1 + d([1 1:nx-1], :, :, 1) - d(:, :, :, 1) + d(:, [1 1:ny-1], :, 2) - d(:,:,:,2);
d(2:end,:,:, 1) = d(2:end,:,:, 1) - d(1:end-1,:,:, 1);
d(1,:,:, 1) = d(1,:,:, 1).*2;
d(:,2:end,:, 2) = d(:,2:end,:, 2) - d(:,1:end-1,:, 2);
d(:,1,:, 2) = d(:,1,:, 2).*2;
u = u - sum(d, 4);

% u1 = u1 - b([1 1:nx-1], :, :, 1) + b(:, :, :, 1) - b(:, [1 1:ny-1], :, 2) + b(:,:,:,2);
b(2:end,:,:, 1) = (b(2:end,:,:, 1) - b(1:end-1,:,:, 1));
b(1,:,:, 1) = b(1,:,:, 1).*2;
b(:,2:end,:, 2) = (b(:,2:end,:, 2) - b(:,1:end-1,:, 2));
b(:,1,:, 2) = b(:,1,:, 2).*2;
u = u + sum(b, 4);

u = u.*(Cl/(1+Cl*4)) + f.*(1/(1+Cl*4));

end