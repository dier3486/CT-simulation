function u = LaplacianSmooth_3Dtanh(img0, lambda, tanhC, u0, Crange, Niter, tol, Zbalance)

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
if nargin<8 || isempty(Zbalance)
    Zbalance = 1.0;
end
imgsize = size(img0);
if length(imgsize) < 3
    % warn! call 2D
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
    delta = zeros(1, Niter, fclass, 'gpuArray');
else
    delta = zeros(1, Niter, fclass);
end

if nargin<4 || isempty(u0)
    u0 = f1;
else
    u0(~s1) = f1(~s1);
end

for ii = 1:Niter
    if ii > 1
        u0 = u;
    end
    u = funG(f1, u0, lambda, tanhC, Zbalance);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
end

if delta(end)>delta(end-1)
    error('Solver failed!');
end

u(~s1) = img0(~s1);

end


function u1 = funG(f, u, lambda, tanhC, Zbalance)
[nx, ny, nz] = size(u);

balance_xy = (3/(2+Zbalance));
balance_z = 3*Zbalance/(2+Zbalance);

du = (u([2:nx nx-1], :, :) + u([2 1:nx-1], :, :) + u(:, [2:ny ny-1], :) + u(:, [2 1:ny-1], :)).*balance_xy ...
     + (u(:, :, [2:nz nz-1]) + u(:, :, [2 1:nz-1])).*balance_z;

% fu = u - tanh((u-f)./tanhC).*tanhC.*(1-tanh((u-f)./tanhC).^2); 
fu = u - tanh((u-f)./tanhC).*tanhC;
% fu = f;
u1 = (du.*lambda + fu)./(1+lambda*6);

end