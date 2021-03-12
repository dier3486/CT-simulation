% bone corr test1

% I know
% [img1, df1, pf1]=CRISrecon('F:\data-Dier.Z\TL\headBC\recon_configure5154.xml', "F:\data-Dier.Z\TL\headBC\2.1611112305154\2.1611112305154.pd");

dataflow = df1;
prmflow = pf1;

Cimage_test = zeros(512,512,2, 'single');
Cimage_test(306:326, 296:316,:) = 1;

gpuDevice;

% FP
Nview = 400;
deltaview = pi/Nview; 
viewangle = linspace(0, pi, Nview+1);
viewangle = single(viewangle(1:Nview));

h_img = prmflow.recon.FOV/prmflow.recon.imagesize;
Nx = prmflow.recon.imagesize;
Ny = prmflow.recon.imagesize;
% Np = 501;
dp = h_img*1.8;
Np = floor(Nx*sqrt(2)*(h_img/dp)/2)*2+1;
d_h = single((-(Np-1)/2:(Np-1)/2)'.*(dp/h_img));

% Cimage1 = reshape(dataflow.image, Nx*Ny, []);
% % Cimage1 = reshape(Cimage_test, Nx*Ny, []);
% Nimg = size(Cimage1, 2);
% Cimage1 = [Cimage1; zeros(1, Nimg)];
% D1 = zeros(Np, Nimg, Nview);
% tic;
% for iview = 1:Nview
%     [inter_a, index_1, index_2, cs_view1] = parallellinearinterp2D(Nx, Ny, d_h, viewangle(iview), Nx*Ny+1);
%     sizesum = size(inter_a);
%     D1(:, :, iview) = squeeze(sum(reshape(Cimage1(index_1(:), :).*(1-inter_a(:)) + ...
%         Cimage1(index_2(:), :).*inter_a(:), [sizesum Nimg]), 2)).*(abs(cs_view1)*h_img);
% end
% toc;

Cimage2 = dataflow.image;
% Cimage2 = Cimage_test;
Nimg = size(Cimage2, 3);
imgindex = repmat(reshape(single(1:Nimg), 1,1,[]), Np, Nx);
Cimage2 = gpuArray(Cimage2);
imgindex = gpuArray(imgindex);
viewangle = gpuArray(viewangle);
d_h = gpuArray(d_h);
Nx = gpuArray(single(Nx));
Ny = gpuArray(single(Ny));
h_img = gpuArray(single(h_img));
viewangle = gpuArray(viewangle);

tic;
D2 = zeros(Np, Nimg, Nview, 'single');
for iview = 1:Nview
    [interpX, interpY, cs_view2] = parallellinearinterp2D2(Nx, Ny, d_h, viewangle(iview));
    Diview = interp3(Cimage2, repmat(interpY,1,1,Nimg), repmat(interpX,1,1,Nimg), imgindex, 'linear', 0); 
    D2(:, :, iview) = gather(squeeze(sum(Diview, 2)).*(abs(cs_view2)*h_img));
end
toc;

D3 = zeros(Np, Nimg, Nview, 'single');
[eta0, zeta0] = meshgrid(d_h);
zetaeta0 = [zeta0(:) eta0(:)];
rze = sqrt(sum(zetaeta0.^2,2));
index_ze = rze<Nx/sqrt(2);
zetaeta0 = zetaeta0(index_ze, :);
Np0 = sum(index_ze);
cosview = cos(viewangle+pi/2);
sinview = sin(viewangle+pi/2);
imgindex = gpuArray(repmat(single(1:Nimg), Np0, 1));
imgcenter = [(Nx+1)/2  (Ny+1)/2];
% 
% tic;
% for iview = 1:Nview/2
%     R = [cosview(iview) sinview(iview); -sinview(iview) cosview(iview)];
%     zetaeta = zetaeta0*R + imgcenter;
%     Diview = zeros(Np^2, Nimg, 'single', 'gpuArray');
%     Diview(index_ze, :) = interp3(Cimage2, repmat(zetaeta(:,2),1,Nimg), repmat(zetaeta(:,1),1,Nimg), imgindex, 'linear', 0);
%     Diview = reshape(Diview, Np, Np, Nimg);
%     D3(:, :, iview) = gather(squeeze(sum(Diview, 2)).*dp);
%     D3(:, :, iview+Nview/2) = gather(squeeze(sum(Diview, 1)).*dp);
% end
% toc;

D4 = zeros(Np, Nimg, Nview, 'single');
[eta1, zeta1] = meshgrid(d_h(1:2:end));
[eta2, zeta2] = meshgrid(d_h(2:2:end));
zetaeta1 = [zeta1(:) eta1(:)];
zetaeta2 = [zeta2(:) eta2(:)];
Np1 = size(zeta1, 1);
Np2 = size(zeta2, 1);
imgindex1 = gpuArray(repmat(single(1:Nimg), Np1^2, 1));
imgindex2 = gpuArray(repmat(single(1:Nimg), Np2^2, 1));

tic;
for iview = 1:Nview/2
    R = [cosview(iview) sinview(iview); -sinview(iview) cosview(iview)];
    % odd
    zetaeta = zetaeta1*R + imgcenter;
    Diview1 = interp3(Cimage2, repmat(zetaeta(:,2),1,Nimg), repmat(zetaeta(:,1),1,Nimg), imgindex1, 'linear', 0);
    Diview1 = reshape(Diview1, Np1, Np1, Nimg);
    D4(1:2:end, :, iview) = gather(squeeze(sum(Diview1, 2)).*(dp*2));
    D4(1:2:end, :, iview+Nview/2) = gather(squeeze(sum(Diview1, 1)).*(dp*2));
    % even
    zetaeta = zetaeta2*R + imgcenter;
    Diview2 = interp3(Cimage2, repmat(zetaeta(:,2),1,Nimg), repmat(zetaeta(:,1),1,Nimg), imgindex2, 'linear', 0);
    Diview2 = reshape(Diview2, Np2, Np2, Nimg);
    D4(2:2:end, :, iview) = gather(squeeze(sum(Diview2, 2)).*(dp*2));
    D4(2:2:end, :, iview+Nview/2) = gather(squeeze(sum(Diview2, 1)).*(dp*2));
end
toc;



