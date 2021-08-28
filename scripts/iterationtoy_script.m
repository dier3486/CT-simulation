% after itertestdata_script.m

% main filter
alpha_filt = 1.5;
Kfilt = gpuArray(filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt));

% FBP
img0a = filterbackproj2D(P0, pbm, Kfilt);
img1a = filterbackproj2D(P1, pbm, Kfilt);

% noise filter
alpha_filt = 2.0;
% Kfilt_BV = filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt);
Kfilt_BV = gpuArray(filterdesign('ram-lak', pbm.Npixel, pbm.delta_d, alpha_filt));

% F^2
Kf2 = Kfilt.*Kfilt_BV.*pi;
% or
Kf1 = Kfilt_BV.*pi;

% Laplace (if use)
% L = splaplace2D(double(pbm.imagesize));

% paramters
mu_BV0 = 0.40;
mu_adap = 0.10;
lambda = 0.03;
mu_F = 2.0;
mu_L = 0.0/mu_F;

% iteration
Niter = 20;
alpha = 0.2;

% vector ini
b1 = img1a;
u = zeros(pbm.imagesize, pbm.imagesize, Niter, 'single');
u(:,:,1) = b1;
Gu = u;
rerr = nan(Niter-1, 1);

% mu ini
u_BV = BregmanTV2D(img1a, mu_BV0, lambda, []);
mu_BV = mu_BV0./(abs(img1a-u_BV).*mu_adap+1);
img1a_tv = u_BV;

% loop
fig = figure;
for ii = 1:Niter-1
    % BV
    u_BV = BregmanTV2D(u(:,:,ii), mu_BV, lambda, u_BV);
    R_BV = (u(:,:,ii) - u_BV);
    mu_BV = mu_BV0./(abs(R_BV).*mu_adap+1);
    Gu(:,:,ii) = u(:,:,ii) - R_BV;
    % style#1 project u;
    % v1 = parallelprojinimageGPU(pbm, u(:,:,ii) + 1i.*R_BV.*mu_F);
    % style#2 project Gu
    v1 = parallelprojinimageGPU(pbm, Gu(:,:,ii) + 1i.*R_BV.*mu_F);
    v = real(v1) + fconv(imag(v1), Kf1);
    r = b1 - filterbackproj2D(v, pbm, Kfilt);
    u(:,:,ii+1) = u(:,:,ii) + gather(r.*alpha);
    
    rerr(ii) = gather(sqrt(sum(r(:).^2)).*(1/pbm.imagesize^2));
    figure(fig);
    plot(1:Niter-1, rerr,'.-');
    axis([0 Niter 0 max(rerr).*1.1]);
    drawnow;
end
Gu(:,:,ii+1) = BregmanTV2D(u(:,:,ii+1), mu_BV, lambda);

function y = fconv(x, K)

sizeX = size(x);
len = size(K, 1);
x(len, 1) = 0;

y = ifft(fft(x).*K);
y = y(1:sizeX,:);

end

function y = LLx(L, x)

xsize = size(x);
y = reshape(L*x(:), xsize);

end

