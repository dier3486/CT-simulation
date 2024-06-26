function u = BregmanTV2DL(img0, mu, mu0, lambda, u0, Crange, Niter, tol)
% Split Bregman method for TV (Laplace version)
% u = BregmanTV2DL(f0, mu, mu0, lambda, u0, Crange, Niter, tol);
% or u = BregmanTV2DL(f0, mu, lambda);
% typicaly, img0 is around 1000 HF, mu=0.1~1, lambda<mu/pi^2


if nargin<6 || isempty(Crange)
    Crange = [-inf inf];
    % Crange = 1000+[-100 100];
end
if nargin<7 || isempty(Niter)
    Niter = 100;
end
if nargin<8 || isempty(tol)
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

% prepare K2
lenxy = cast(2.^(floor(log2(imgsize(1:2)))+1), fclass);
kx = [0:lenxy(1)/2  -lenxy(1)/2+1:-1].*(pi*2/lenxy(1));
ky = [0:lenxy(2)/2  -lenxy(2)/2+1:-1].*(pi*2/lenxy(2));
[Ky, Kx] = meshgrid(ky, kx);
K2 = 1./(mu0-(Ky.^2 + Kx.^2).*lambda);

% gpuArray
if GPUonoff
    b = zeros([imgsize 4], fclass, 'gpuArray');
    d = zeros([imgsize 4], fclass, 'gpuArray');
    delta = zeros(1, Niter, fclass, 'gpuArray');
    K2 = gpuArray(K2);
else
    b = zeros([imgsize 4], fclass);
    d = zeros([imgsize 4], fclass);
    delta = zeros(1, Niter, fclass);
end

if nargin<4 || isempty(u0)
    u0 = f1;
else
    u0(~s1) = f1(~s1);
    [d, b] = fundbyub(u0, b, lambda);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, b, d, mu, mu0, lambda, K2, lenxy);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(u, b, lambda);
end

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

% [nx, ny, nz] = size(u);

% du = zeros(nx, ny, nz, 2);
du = b0.*0;
du(1:end-1,:,:,1) = u(2:end,:,:) - u(1:end-1, :,:);
du(:,1:end-1,:,2) = u(:, 2:end,:) - u(:, 1:end-1,:);

du(1:end-1,1:end-1,:,3) = u(2:end,2:end,:) - u(1:end-1, 1:end-1,:);
du(2:end,1:end-1,:,4) = u(1:end-1,2:end,:) - u(2:end, 1:end-1,:);


b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);
d1 = max(absb0-1/lambda, 0).*(b0./absb0);
% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, mu0, lambda, K2, lenxy)
[nx, ny, ~] = size(u);

% u1 = u([2:nx nx], :, :) + u([1 1:nx-1], :, :) + u(:, [2:ny ny], :) + u(:, [1 1:ny-1], :);
% u1 = u1 + d([1 1:nx-1], :, :, 1) - d(:, :, :, 1) + d(:, [1 1:ny-1], :, 2) - d(:,:,:,2);
% u1 = u1 - b([1 1:nx-1], :, :, 1) + b(:, :, :, 1) - b(:, [1 1:ny-1], :, 2) + b(:,:,:,2);
% u1 = u1.*(lambda./(mu+4*lambda)) + f.*(mu./(mu+4*lambda));

u1 = f.*mu - (mu-mu0 ).*u;
u1 = u1 + (d([1 1:nx-1], :, :, 1) - d(:, :, :, 1) + d(:, [1 1:ny-1], :, 2) - d(:,:,:,2)).*lambda;
u1 = u1 + (-b([1 1:nx-1], :, :, 1) + b(:, :, :, 1) - b(:, [1 1:ny-1], :, 2) + b(:,:,:,2)).*lambda;

% u1 = u1 + (d([1 1:nx-1], :, :, 3) - d(:, :, :, 3) + d(:, [1 1:ny-1], :, 4) - d(:,:,:,4)).*lambda./2;
% u1 = u1 + (-b([1 1:nx-1], :, :, 3) + b(:, :, :, 3) - b(:, [1 1:ny-1], :, 4) + b(:,:,:,4)).*lambda./2;

u2 = ifft2(fft2(u1, lenxy(1), lenxy(2)).*K2, 'symmetric');

u1 = u2(1:nx, 1:ny, :);

end