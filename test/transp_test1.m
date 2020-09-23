[N, ~, Nslice] = size(Cimage);
h = 1;
hond = h/pb.delta_d;
centerond = [0 0];
theta = pb.viewangle;
Np = pb.Np;
Nview = length(theta);
ctrIdx = (Np+1)/2;

% Define the x & y axes for the reconstructed image
[x, y] = ndgrid(-(N-1)/2 : (N-1)/2);
x = x(:).*hond;
y = y(:).*hond;

% Generate trignometric tables
costheta = cos(theta);
sintheta = sin(theta);

p = zeros(Np, Nslice, Nview);
p1 = zeros(Np, Nslice, Nview);
sizep = size(p);

% Allocate memory for the image
img = zeros(N,N, Nslice, 'like',p);
Cimage = reshape(Cimage, N*N, Nslice);

% interp
for iview=1:length(theta)
    % fill 0 and nan
%     proj = [p(:, :, iview); zeros(1, Nslice); nan(1, Nslice)];
    proj = zeros(Np+2, Nslice);
    proj1 = zeros(Np+2, Nslice);
    % projection sample
    t = (-x.*sintheta(iview) + y.*costheta(iview)) + ctrIdx;
    % interpolation index and alpha
    tindex = floor(t);
    t_alpha = t - tindex;
    % index+1
    tindex_p1 = tindex + 1;
    % edge
    tindex(tindex==0) = Np+1;
    tindex(tindex<0) = Np+2;
    tindex(tindex>Np) = Np+1;
    tindex_p1(tindex_p1<=0) = Np+2;
    tindex_p1(tindex_p1>Np+1) = Np+2;
%     % slice shift
%     tindex = tindex + (Np+2).*(0 : Nslice-1);
%     tindex_p1 = tindex_p1 + (Np+2).*(0 : Nslice-1);
%     % get sample value by interpolation
%     projContrib = proj(tindex(:)).*(1-t_alpha(:)) + proj(tindex_p1(:)).*t_alpha(:);
%     % add to image
%     img = img + reshape(projContrib, N, N, Nslice);
    % 1
    for ii = 1:N^2
        proj(tindex(ii), :) = proj(tindex(ii), :) +  Cimage(ii,:).*(1-t_alpha(ii));
        proj(tindex_p1(ii), :) = proj(tindex_p1(ii), :) +  Cimage(ii,:).*t_alpha(ii);
        proj1(tindex(ii), :) = proj1(tindex(ii), :) +  (1-t_alpha(ii));
        proj1(tindex_p1(ii), :) = proj1(tindex_p1(ii), :) +  t_alpha(ii);
    end
    p(:, :, iview) = proj(1:Np, :);
    p1(:, :, iview) = proj1(1:Np, :);
end

Cimage = reshape(Cimage, N, N, Nslice);
% img = img.*(pi/length(theta)/2);
