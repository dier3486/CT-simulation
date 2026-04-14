function [u, b, d] = BregmanTV3D_tanh(img0, mu, lambda, u0, Crange, Niter, tol, Zbalance, tanhSigma, tanhC)
% Split Bregman method for 3D TV
% u = BregmanTV3D(f0, mu, lambda, u0, Crange, Niter, tol);
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
if nargin<8 || isempty(Zbalance)
    Zbalance = 1.0;
end
if nargin<9 || isempty(tanhSigma)
    tanhSigma = 1e-2;
    tanhC = 2;
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
    b = zeros([imgsize 3], fclass, 'gpuArray');
    d = zeros([imgsize 3], fclass, 'gpuArray');
    delta = zeros(1, Niter, fclass, 'gpuArray');
else
    b = zeros([imgsize 3], fclass);
    d = zeros([imgsize 3], fclass);
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
    u = funG(f1, u0, b, d, mu, lambda, Zbalance, tanhSigma, tanhC);
    delta(ii) = sqrt(sum(abs(u(:)-u0(:)).^2)./N1);
    if delta(ii)<tol
        break;
    end
    [d, b] = fundbyub(u, b, lambda);
%     disp(delta(ii));
end

disp([ii delta(ii)]);

if Niter>1 && delta(end)>delta(end-1)
    error('TV failed!');
end

% u = u./1000;
u(~s1) = img0(~s1);
% delta(ii)

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
d1 = max(absb0-1/lambda, 0).*(b0./absb0);

% d1 = fillmissing(d1, 'constant', 0);
% d1(isnan(d1)) = 0;

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, lambda, Zbalance, tanhSigma, tanhC)
[nx, ny, nz] = size(u);

balance_xy = (3/(2+Zbalance));
balance_z = 3/(2/Zbalance+1);
balance = reshape([balance_xy balance_xy balance_z], 1, 1, 1, 3);

% lambda_u = u.*0 + lambda;
% lambda_u(2:nx-1, 2:ny-1, 2:nz-1) = lambda_u(2:nx-1, 2:ny-1, 2:nz-1).*1.5;

% u1 = u([2:nx nx], :, :) + u([1 1:nx-1], :, :) + u(:, [2:ny ny], :) + u(:, [1 1:ny-1], :) ...
%      + u(:, :, [2:nz nz]) + u(:, :, [1 1:nz-1]);
u1 = (u([2:nx nx-1], :, :) + u([2 1:nx-1], :, :) + u(:, [2:ny ny-1], :) + u(:, [2 1:ny-1], :)).*balance_xy ...
     + (u(:, :, [2:nz nz-1]) + u(:, :, [2 1:nz-1])).*balance_z;

% u1 = u1 + d([1 1:nx-1], :, :, 1) - d(:, :, :, 1) + d(:, [1 1:ny-1], :, 2) - d(:, :, :, 2) ...
%      + d(:, :, [1 1:nz-1], 3) - d(:, :, :, 3);
d(2:end,:,:, 1) = d(2:end,:,:, 1) - d(1:end-1,:,:, 1);
d(1,:,:, 1) = d(1,:,:, 1).*2;
d(:,2:end,:, 2) = d(:,2:end,:, 2) - d(:,1:end-1,:, 2);
d(:,1,:, 2) = d(:,1,:, 2).*2;
d(:,:,2:end, 3) = d(:,:,2:end, 3) - d(:,:,1:end-1, 3);
d(:,:,1, 3) = d(:,:,1, 3).*2;
d = d.*balance;
u1 = u1 - sum(d, 4);

% u1 = u1 - b([1 1:nx-1], :, :, 1) + b(:, :, :, 1) - b(:, [1 1:ny-1], :, 2) + b(:,:,:,2) ...
%      - b(:, :, [1 1:nz-1], 3) + b(:, :, :, 3);
b(2:end,:,:, 1) = (b(2:end,:,:, 1) - b(1:end-1,:,:, 1));
b(1,:,:, 1) = b(1,:,:, 1).*2;
b(:,2:end,:, 2) = (b(:,2:end,:, 2) - b(:,1:end-1,:, 2));
b(:,1,:, 2) = b(:,1,:, 2).*2;
b(:,:,2:end, 3) = (b(:,:,2:end, 3) - b(:,:,1:end-1, 3));
b(:,:,1, 3) = b(:,:,1, 3).*2;
b = b.*balance;
u1 = u1 + sum(b, 4);

f1 = hardentanh2(f-u, tanhSigma, tanhC) + u;


u2 = u1.*(lambda./(mu+lambda*6)) + f1.*(mu./(mu+lambda*6));
% u1 = u1.*(lambda./(mu+lambda*8)) + f.*(mu./(mu+lambda*8));

% f2 = f-u;
% f2(abs(f2)<tanhC) = tanhC.*sign(f2(abs(f2)<tanhC));
% f2 = f2+u;
% 
% 1;
 
% v = u;
% alpha = 1.0;
% N = 10;
% dr = zeros(N, nz);
% for ii = 1:N
%     r = (hardentanh2(f-v, tanhSigma, tanhC)+v).*(mu./(mu+lambda*6)) + u1.*(lambda./(mu+lambda*6)) - v;
%     v = v + r.*alpha;
%     dr(ii, :) = sqrt(sum(reshape(r.^2, [], nz)))./nx./ny;
% end
% 
% beta = 0.0;
% 
% u1 = u2 + (v-u2).*beta;

% u1 = u1.*(lambda./(mu+lambda*6)) + f2.*(mu./(mu+lambda*6));

u1 = u2;

end