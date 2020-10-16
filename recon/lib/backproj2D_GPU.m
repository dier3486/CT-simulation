function img = backproj2D_GPU(p, theta, ctrIdx, hond, N, centerond, maxR)
% back projection 
% img = backproj2D_GPU(p, theta, ctrIdx, hond, N, centerond, maxR)
% like backproj2D_2

if nargin<6
    centerond = [0 0];
end
if nargin<7
    maxR = inf;
end

pclass = class(p);

% size of the p, Np*Nslice*Nview
sizep = size(p);
Np = sizep(1);
if length(sizep) > 2
    Nslice = sizep(2);
else
    Nslice = 1;
    p = reshape(p, sizep(1), 1, sizep(2));
end

% Define the x & y axes for the reconstructed image
[x, y] = ndgrid(-(N-1)/2 : (N-1)/2);
x = x(:).*hond - centerond(:, 1)';
y = y(:).*hond - centerond(:, 2)';
if size(x, 2)==1 && Nslice>1
    x = repmat(x, 1, Nslice);
    y = repmat(y, 1, Nslice);
end
Sxy = x(:).^2 + y(:).^2 <= maxR.^2;
x = x(Sxy);
y = y(Sxy);
Nxy = sum(Sxy);
% slice index shift
sliceshift = repmat(0 : Nslice-1, N^2, 1).*Np;
sliceshift = sliceshift(Sxy);

% Generate trignometric tables
costheta = cos(theta(:)');
sintheta = sin(theta(:)');

% to gpu
x = gpuArray(cast(x, pclass));
y = gpuArray(cast(y, pclass));
costheta = gpuArray(cast(costheta, pclass));
sintheta = gpuArray(cast(sintheta, pclass));
ctrIdx = gpuArray(cast(ctrIdx, pclass));
% sliceshift = gpuArray(cast(sliceshift, pclass));
Np = gpuArray(cast(Np, pclass));
% Nslice = gpuArray(cast(Nslice, pclass));
% N = gpuArray(cast(N, pclass));
% proj = gpuArray([zeros(Np+1, Nslice, pclass); nan(1, Nslice, pclass)]);
% proj = gpuArray(zeros(Np, Nslice, pclass));
% Allocate memory for the image
img_GPU = gpuArray(zeros(Nxy, 1, pclass));

Nview = length(theta);
maxview = 4;
% warn: too large maxview could: 1. out of GPU memory 2. out of single accuracy on the index,
viewloop = ceil(Nview/maxview);
proj = gpuArray(zeros(Np*Nslice, maxview, pclass));
sliceshift = gpuArray(cast(sliceshift(:) + (0:maxview-1).*(Np*Nslice), pclass));
for iv = 1:viewloop
    % views
    v1 = (iv-1)*maxview+1;
    if iv<viewloop
        v2 = v1+maxview-1;
        Nv = maxview;
    else
        v2 = Nview;
        Nv = v2-v1+1;
    end
    % p -> GPU
    proj(:, 1:Nv) = reshape(p(:, :, v1:v2), Np*Nslice, Nv);
    % projection sample
    t = (-x.*sintheta(v1:v2) + y.*costheta(v1:v2)) + ctrIdx;
     % interpolation index and alpha
    tindex = floor(t);
    t_alpha = t - tindex;
    % slice shift
    tindex = tindex + sliceshift;
    % index+1
%     tindex_p1 = tindex + 1;
    % get sample value by interpolation
    projContrib = proj(tindex).*(1-t_alpha) + proj(tindex+1).*t_alpha;
    % add to image
    img_GPU = img_GPU + sum(projContrib, 2);
end

img = zeros(N, N, Nslice, pclass);
img(Sxy) = gather(img_GPU).*(pi/length(theta)/2);

end