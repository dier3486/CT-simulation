function [u, b] = BregmanTV3D_helicaliter(u, mu, Cl, b0, Crange, Niter, tol, Zbalance)
% Split Bregman method for 3D TV, tmp version for helical iteration recon
% [u, b] = BregmanTV3D_helicaliter(u, mu, Cl, b0, Crange, Niter, tol, Zbalance)
% or u = BregmanTV3D_axialiter(u, mu, Cl);
% where u was u + iGu
% the Cl is lambda./mu
% typicaly, img0 is around 1000 HF, mu=0.1~1, lambda<mu/4

if nargin<4
    b0 = [];
end
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
imgsize = size(u);
if length(imgsize) < 3
    % warn! call 2D
    imgsize = [imgsize 1];
end

GPUonoff = isa(u, 'gpuArray');
fclass = classGPU(u);

f1 = real(u);
s1 = (f1>=Crange(1)) & (f1<=Crange(2));
% N1 = sum(s1(:));
% f1(f1<Crange(1)) = Crange(1);
% f1(f1>Crange(2)) = Crange(2);
N1 = size(u(:), 1);

u0 = imag(u);
if ~isempty(b0)
    [d, b] = fundbyub(u0, b0, mu.*Cl);
else
    if GPUonoff
        b = zeros([imgsize 3], fclass, 'gpuArray');
        d = zeros([imgsize 3], fclass, 'gpuArray');
    else
        b = zeros([imgsize 3], fclass);
        d = zeros([imgsize 3], fclass);
    end
end

if GPUonoff
    delta = zeros(1, Niter, fclass, 'gpuArray');
else
    delta = zeros(1, Niter, fclass);
end

for ii = 1:Niter
    if ii > 1
        u0 = imag(u);
    end
    u = real(u) + 1.0i.*funG(real(u), u0, b, d, Cl, Zbalance);
    delta(ii) = sqrt(sum((imag(u(:))-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(imag(u), b, mu.*Cl);
end

if delta(end)>delta(end-1)
    error('TV failed!');
end

u(~s1) = f1(~s1);
% delta(ii)
% disp([ii delta(ii)]);
end


function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)
% d_y^{k+1} = shrink(D_yu^{k+1}+b_y^{k}, 1/\lambda)
% d_z^{k+1} = shrink(D_zu^{k+1}+b_z^{k}, 1/\lambda)

[nx, ny, nz] = size(u);
% du = zeros(nx, ny, nz, 3);

du = b0.*0;
% du(1:end-1,:,:, 1) = u(2:end, :,:) - u(1:end-1, :,:);
% du(:,1:end-1,:, 2) = u(:, 2:end,:) - u(:, 1:end-1,:);
% du(:,:,1:end-1, 3) = u(:,:, 2:end) - u(:,:, 1:end-1); 

du(:,:,:, 1) = u([2:nx  nx], :,:) - u([1:nx-1  nx-1], :,:);
du(:,:,:, 2) = u(:, [2:ny  ny],:) - u(:, [1:ny-1  ny-1],:);
du(:,:,:, 3) = (u(:,:, [2:nz  nz]) - u(:,:, [1:nz-1  nz-1]));

du(end,:,:, 1) = -du(end,:,:, 1);
du(:,end,:, 2) = -du(:,end,:, 2);
du(:,:,end, 3) = -du(:,:,end, 3);

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);   % to replace fillmissing
d1 = max(absb0-1./lambda, 0).*(b0./absb0);

% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function u = funG(f, u, b, d, Cl, Zbalance)
[nx, ny, nz] = size(u);

balance_xy = 3/(2+Zbalance);
balance_z = 3/(2/Zbalance+1);
balance = reshape([balance_xy balance_xy balance_z], 1, 1, 1, 3);

% u1 = u([2:nx nx], :, :) + u([1 1:nx-1], :, :) + u(:, [2:ny ny], :) + u(:, [1 1:ny-1], :) ...
%      + u(:, :, [2:nz nz]) + u(:, :, [1 1:nz-1]);
u(:,:, 2:end) = (u([2:nx nx-1], :, 2:end) + u([2 1:nx-1], :, 2:end) + ...
    u(:, [2:ny ny-1], 2:end) + u(:, [2 1:ny-1], 2:end)).*balance_xy + ...
    (u(:, :, [3:nz nz-1]) + u(:, :, [1:nz-1])).*balance_z;

% u1 = u1 + d([1 1:nx-1], :, :, 1) - d(:, :, :, 1) + d(:, [1 1:ny-1], :, 2) - d(:, :, :, 2) ...
%      + d(:, :, [1 1:nz-1], 3) - d(:, :, :, 3);
d(2:end,:,:, 1) = d(2:end,:,:, 1) - d(1:end-1,:,:, 1);
d(1,:,:, 1) = d(1,:,:, 1).*2;
d(:,2:end,:, 2) = d(:,2:end,:, 2) - d(:,1:end-1,:, 2);
d(:,1,:, 2) = d(:,1,:, 2).*2;
d(:,:,2:end, 3) = d(:,:,2:end, 3) - d(:,:,1:end-1, 3);
% d(:,:,1, 3) = d(:,:,1, 3).*2;
d = d.*balance;
u(:,:, 2:end) = u(:,:, 2:end) - sum(d(:,:, 2:end, :), 4);

% u1 = u1 - b([1 1:nx-1], :, :, 1) + b(:, :, :, 1) - b(:, [1 1:ny-1], :, 2) + b(:,:,:,2) ...
%      - b(:, :, [1 1:nz-1], 3) + b(:, :, :, 3);
b(2:end,:,:, 1) = (b(2:end,:,:, 1) - b(1:end-1,:,:, 1));
b(1,:,:, 1) = b(1,:,:, 1).*2;
b(:,2:end,:, 2) = (b(:,2:end,:, 2) - b(:,1:end-1,:, 2));
b(:,1,:, 2) = b(:,1,:, 2).*2;
b(:,:,2:end, 3) = (b(:,:,2:end, 3) - b(:,:,1:end-1, 3));
% b(:,:,1, 3) = b(:,:,1, 3).*2;
b = b.*balance;
u(:,:, 2:end) = u(:,:, 2:end) + sum(b(:,:, 2:end, :), 4);
 
u(:,:, 2:end) = u(:,:, 2:end).*(Cl/(1+Cl*6)) + f(:,:, 2:end).*(1/(1+Cl*6));

end