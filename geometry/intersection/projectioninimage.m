function D = projectioninimage(focalspot, detectposition, Cimage, viewangle, couch, method)
% 3D projection on image(s)
% D = projectioninimage(focalspot, detectposition, Cimage, viewangle, couch, method)
% shall be faster than the genernal function intersection in as
% intersection.m

if nargin<4
    viewangle = 0;
end

if nargin<5
    couch = zeros(size(viewangle));
elseif isempty(couch)
    couch = zeros(size(viewangle));
end

if nargin<6
    method = 'lineinsection';
elseif isempty(method)
    method = 'lineinsection';
end

% image(s) grid
[Nx, Ny, Nz] = size(Cimage);
h = 1.0;
Xgrid = (-Nx/2:Nx/2).*h;
Ygrid = -Ny/2:Ny/2;
Zgrid = -Nz/2:Nz/2;
% N
Ndet = size(detectposition, 1);
Nview = size(viewangle(:), 1);
% Lxy is the length of AB on xy plane
Lxy = sqrt((detectposition(:,1)-focalspot(1)).^2 + (detectposition(:,2)-focalspot(2)).^2);
% Lsec = L/Lxy
Lsec = sqrt(sum((detectposition - repmat(focalspot, Ndet, 1)).^2, 2))./Lxy;
% d is the distance from ISO to AB
d = (focalspot(1).*detectposition(:,2)-detectposition(:,1).*focalspot(2))./Lxy;
Lmid = sqrt(focalspot(1).^2+focalspot(2).^2-d.^2);
% Zctg is the ctg(theta_z) = Lxy/Z_AB
Zctg = Lxy./(detectposition(:,3)-focalspot(3));
% Z_A is A(:,3);
Z_A = focalspot(3);
% angles
theta0 = atan2(detectposition(:,2)-focalspot(2), detectposition(:,1)-focalspot(1));

switch method
    case 'lineinsection'
        D = zeros(Ndet, Nview);
        for iview = 1:Nview
            theta = theta0+viewangle(iview);
            Z_Ai = Z_A - couch(iview);
            [dt, Vindex] = linesinimage3D(theta, d, [], Lmid, Z_Ai, Zctg, Xgrid, Ygrid, Zgrid);
            Cimage = [Cimage(:); 0];
            D(:, iview) = sum(dt.*Cimage(Vindex), 2).*Lsec;
        end
    otherwise
        D = [];
end

return