% after itertestdata_script.m

% part a
pbm_a = pbm;
pbm_a.Nview = 500;
pbm_a.viewangle = pbm.viewangle(1:500);
P1a = P1(:,1:500);

% part b
Npart = 900;
pbm_b = pbm;
pbm_b.Nview = Npart;
pbm_b.viewangle = pbm.viewangle(1:Npart);
P1b = P1(:,1:Npart);

% filters
Klr = gpuArray(filterdesign('ram-lak', pbm.Npixel, pbm.delta_d, 2.0));

% main filter
alpha_lr = 0;
alpha_filt = 1.5;
Kfilt = gpuArray(filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt));
Kfilt = Kfilt.*(1-alpha_lr) + Klr.*alpha_lr;

% FBP
img1a = filterbackproj2D(P1a, pbm_a, Kfilt);
img1b = filterbackproj2D(P1b, pbm_b, Kfilt);

% noise filter
alpha_filt = 1.5;
alpha_lr = 1.0;
Kfilt_BV = filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt);
% Kfilt_BV = gpuArray(filterdesign('ram-lak', pbm.Npixel, pbm.delta_d, alpha_filt));
Kfilt_BV = Kfilt_BV.*(1-alpha_lr) + Klr.*alpha_lr;

% F^2
Kf2 = Kfilt.*Kfilt_BV.*pi;
% or
Kf1 = Kfilt_BV.*pi;

% Laplace (if use)
L = splaplace2D(double(pbm.imagesize));

% D weight
sigma_ERF = 0.5;
sigma1 = 1./sqrt(exp(-P1b.*1e-3.*mu_e).*Ndose)./mu_e.*1e3./sqrt(pbm.Nview*pbm.imagesize);
D = gpuArray(ones(size(P1), classGPU(P1b)));
D(:, 1:Npart) = erf(sigma_ERF./sigma1);
D(D<0.05) = 0.05;

% paramters
mu_BV0 = 0.30;
mu_adap = 0.15;
lambda = 0.01;
mu_F = 1.0;
mu_L = 0.03/mu_F;

% iteration
Niter = 20;
alpha = 0.2;

% vector ini
b1 = img1b;
u = zeros(pbm.imagesize, pbm.imagesize, Niter, 'single');
u(:,:,1) = img1a;
Gu = u;
rerr = nan(Niter-1, 1);

% mu ini
u_BV = BregmanTV2D(img1a, mu_BV0, lambda, []);
mu_BV = mu_BV0./(abs(img1a-u_BV).*mu_adap+1);
% mu_BV = mu_BV0.*(1-erf(abs(img1a-u_BV).*mu_adap));   mu_BV(mu_BV<0.05) = 0.05;
% to compare 
img1a_tv0 = u_BV;
img1a_tv1 = BregmanTV2D(img1a, mu_BV, lambda, u_BV);

% loop
fig = figure;
for ii = 1:Niter-1
    % BV
    u_BV = BregmanTV2D(u(:,:,ii), mu_BV, lambda, u_BV);
    R_BV = (u(:,:,ii) - u_BV);
    mu_BV = mu_BV0./(abs(R_BV).*mu_adap+1);
%     mu_BV = mu_BV0.*(1-erf(abs(img1a-u_BV).*mu_adap));   mu_BV(mu_BV<0.05) = 0.05;
    Gu(:,:,ii) = u(:,:,ii) - R_BV;
    % L
    if mu_L>0
        R_BV = R_BV + LLx(L, Gu(:,:,ii)).*mu_L;
%         R_BV = R_BV + LLx(L, (Gu(:,:,ii)+u(:,:,ii))./2).*mu_L;
    end
    
    % style#1 project u;
    % v1 = parallelprojinimageGPU(pbm, u(:,:,ii));
    % style#2 project Gu
%     v1 = parallelprojinimageGPU(pbm, Gu(:,:,ii));
    % style#3 conjugate
    v1 = parallelprojinimageGPU(pbm_b, (Gu(:,:,ii)+u(:,:,ii))./2);
    % v2
    v2 = parallelprojinimageGPU(pbm, R_BV.*mu_F)./D;
    % FBP
    r1 = filterbackproj2D(v1, pbm_b, Kfilt);
    r2 = filterbackproj2D(v2, pbm, Kf2);
    
%     v = real(v1) + fconv(imag(v1), Kf1)./D;
    r = b1 - r1 - r2;
    u(:,:,ii+1) = u(:,:,ii) + gather(r.*alpha);
    
    rerr(ii) = gather(sqrt(sum(r(:).^2)).*(1/pbm.imagesize^2));
    figure(fig);
    plot(1:Niter-1, rerr,'.-');
    axis([0 Niter 0 max(rerr).*1.1]);
    drawnow;
end
Gu(:,:,ii+1) = BregmanTV2D(u(:,:,ii+1), mu_BV, lambda, u_BV);

function y = fconv(x, K)

sizeX = size(x);
len = size(K, 1);
x(len, 1) = 0;

y = ifft(fft(x).*K);
y = y(1:sizeX,:);

end

function y = LLx(L, x)

xclass = classGPU(x);
xsize = size(x);
y = cast(reshape(L*double(x(:)), xsize), xclass);

end

