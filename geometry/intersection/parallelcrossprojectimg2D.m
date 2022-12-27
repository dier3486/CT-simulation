function D = parallelcrossprojectimg2D(Cimage, d_h, zetaangle, imgcenter)
% D = parallelcrossprojectimg2D(Cimage, d./h, zetaangle, imgcenter).*delta_d;
% I know zetaangle = viewangle+pi/2
% WARN: known bug: the Cimage shall be Cimage.'

if nargin<4
    imgcenter = [0 0];
end

% GPU?
gpuonoff = isa(Cimage,'gpuArray');
% becuase it will call interp2 or interp3, strongly suggest to employ GPU

dataclass = classGPU(Cimage);

Nview = size(zetaangle(:), 1);
% I know the Nview is half of the view number which covered pi/2
Np = size(d_h(:), 1);
[Nx, Ny, Nimg] = size(Cimage);

% zeta-eta grid, odd and even
[eta1, zeta1] = meshgrid(d_h(1:2:end));
[eta2, zeta2] = meshgrid(d_h(2:2:end));
zetaeta1 = [zeta1(:) eta1(:)];
zetaeta2 = [zeta2(:) eta2(:)];
Np1 = size(zeta1, 1);
Np2 = size(zeta2, 1);
% image center and image index for interpolation
if gpuonoff
    imgindex1 = gpuArray(repmat(cast(1:Nimg, dataclass), Np1^2, 1));
    imgindex2 = gpuArray(repmat(cast(1:Nimg, dataclass), Np2^2, 1));
    imgcenter = gpuArray(cast(imgcenter + [(Nx+1)/2  (Ny+1)/2], dataclass));
else
    imgindex1 = repmat(cast(1:Nimg, dataclass), Np1^2, 1);
    imgindex2 = repmat(cast(1:Nimg, dataclass), Np2^2, 1);
    imgcenter = cast(imgcenter + [(Nx+1)/2  (Ny+1)/2], dataclass);
end
% cos sin
cosview = cos(zetaangle);
sinview = sin(zetaangle);

% ini D
D = zeros(Np, Nimg, Nview*2, classGPU(Cimage));
% loop view
for iview = 1:Nview
    % rotation matrix
    R = [cosview(iview) sinview(iview); -sinview(iview) cosview(iview)];
    % odd
    zetaeta = zetaeta1*R + imgcenter;
    Diview1 = interp3(Cimage, repmat(zetaeta(:,2),1,Nimg), repmat(zetaeta(:,1),1,Nimg), imgindex1, 'linear', 0);
    Diview1 = reshape(Diview1, Np1, Np1, Nimg);
    D(1:2:end, :, iview) = gather(squeeze(sum(Diview1, 2)).*2);
    D(1:2:end, :, iview+Nview) = gather(squeeze(sum(Diview1, 1)).*2);
    % even
    zetaeta = zetaeta2*R + imgcenter;
    Diview2 = interp3(Cimage, repmat(zetaeta(:,2),1,Nimg), repmat(zetaeta(:,1),1,Nimg), imgindex2, 'linear', 0);
    Diview2 = reshape(Diview2, Np2, Np2, Nimg);
    D(2:2:end, :, iview) = gather(squeeze(sum(Diview2, 2)).*2);
    D(2:2:end, :, iview+Nview) = gather(squeeze(sum(Diview2, 1)).*2);
end

end