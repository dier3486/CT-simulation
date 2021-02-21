function [v, w] = focalresample(focalpos, focalsize, Nresmp)
% detector resample

Nx = Nresmp(1);
if size(Nresmp(:),1)>1
    Nz = Nresmp(2);
else
    Nz = 1;
end
Nfocal = size(focalpos, 1);

% sample nodes and weight
[px, wx] = lgwt(Nx, -0.5, 0.5);
[pz, wz] = lgwt(Nz, -0.5, 0.5);
w = wx*wz';
w = w(:);
Np = size(w, 1);
[PZ, PX] = meshgrid(pz,px);
% Vx Vz
Vx = repmat([focalsize(1) 0 0], Nfocal, 1);
Vz = repmat([0 0 focalsize(2)], Nfocal, 1);

% v
v = focalpos(:) + (Vx(:)*PX(:)' + Vz(:)*PZ(:)');
v = reshape(v, [], 3*Np);

end