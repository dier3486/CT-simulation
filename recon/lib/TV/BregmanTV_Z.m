function u = BregmanTV_Z(img0, mu, lambda, u0, Crange, Niter, tol)
% Split Bregman method for 3D-Z TV
% u = BregmanTV_Z(f0, mu, lambda, u0, Crange, Niter, tol);
% or u = BregmanTV3D(f0, mu, lambda);
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
    error('2D image can not besmoothed on z-direction');
end

f1 = img0;
s1 = (f1>=Crange(1)) & (f1<=Crange(2));
N1 = sum(s1(:));
f1(f1<Crange(1)) = Crange(1);
f1(f1>Crange(2)) = Crange(2);

GPUonoff = isa(img0, 'gpuArray');
fclass = classGPU(img0);

if GPUonoff
    b = zeros(imgsize, fclass, 'gpuArray');
    d = zeros(imgsize, fclass, 'gpuArray');
    delta = zeros(1, Niter, fclass, 'gpuArray');
else
    b = zeros(imgsize, fclass);
    d = zeros(imgsize, fclass);
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
    u = funG(f1, u0, b, d, mu, lambda);
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
% d_z^{k+1} = shrink(D_zu^{k+1}+b_z^{k}, 1/\lambda)

nz = size(u, 3);

du = (u(:,:, [2:nz  nz]) - u(:,:, [1:nz-1  nz-1]));
du(:,:,end) = -du(:,:,end);

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);   % to replace fillmissing
d1 = max(absb0-1/lambda, 0).*(b0./absb0);

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, lambda)
nz = size(u, 3);

u1 = u(:, :, [2:nz nz-1]) + u(:, :, [2 1:nz-1]);

d(:,:,2:end) = d(:,:,2:end) - d(:,:,1:end-1);
d(:,:,1) = d(:,:,1).*2;
u1 = u1 - d;

b(:,:,2:end) = b(:,:,2:end) - b(:,:,1:end-1);
b(:,:,1) = b(:,:,1).*2;
u1 = u1 + b;

f1 = hardentanh2(f-u, 1e-2, 2) + u;

u1 = u1.*(lambda./(mu+lambda*2)) + f1.*(mu./(mu+lambda*2));

end