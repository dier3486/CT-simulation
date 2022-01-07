function [x,y,z,Sxy,costheta,sintheta] = backproj2Dprepare(N, Nslice, hond, centerond, maxR, theta, pclass)

% Define the x & y axes for the reconstructed image
[x, y] = ndgrid(-(N-1)/2 : (N-1)/2);
x = x(:).*hond - centerond(:, 1)';
y = y(:).*hond - centerond(:, 2)';
if size(x, 2)==1 && Nslice>1
    x = repmat(x, 1, Nslice);
    y = repmat(y, 1, Nslice);
end
Sxy = any(x.^2 + y.^2 <= maxR.^2, 2);

x = x(Sxy, :);
y = y(Sxy, :);
Nxy = sum(Sxy);
% z (slice)
z = repmat(1:Nslice, Nxy, 1);

% Generate trignometric tables
costheta = cos(theta(:)');
sintheta = sin(theta(:)');

% to gpu
x = gpuArray(cast(x, pclass));
y = gpuArray(cast(y, pclass));
z = gpuArray(cast(z, pclass));
costheta = gpuArray(cast(costheta, pclass));
sintheta = gpuArray(cast(sintheta, pclass));
% ctrIdx = gpuArray(cast(ctrIdx, pclass));