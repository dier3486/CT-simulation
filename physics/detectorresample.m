function [v, w] = detectorresample(detector, Nresmp)
% detector resample

Nx = Nresmp(1);
if size(Nresmp(:),1)>1
    Nz = Nresmp(2);
else
    Nz = 1;
end

% sample nodes and weight
[px, wx] = lgwt(Nx, -0.5, 0.5);
[pz, wz] = lgwt(Nz, -0.5, 0.5);
w = wx*wz';
w = w(:);
Np = size(w, 1);
[PZ, PX] = meshgrid(pz,px);

% x
Vx = zeros(size(detector.normvector));
Vx(:, 1) = -detector.normvector(:, 2);
Vx(:, 2) = detector.normvector(:, 1);
Vx = normr(Vx);
% z
Vz = cross(detector.normvector, Vx);
% length scale
Vx = Vx.*detector.edgelength(:, 1);
Vz = Vz.*detector.edgelength(:, 2);

% v
v = detector.position(:) + (Vx(:)*PX(:)' + Vz(:)*PZ(:)');
v = reshape(v, [], 3*Np);

end