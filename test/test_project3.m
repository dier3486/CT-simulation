% test script for ideal detector position
% on grid image

addpath(genpath('../'));

addpath(genpath('../'));

% focalpos_ideal = [0, -550, 0];
% path AB (as A is the focal spot, B is the detector pixels)
B = [100, -500, 0.1; 100, -500, 10];
A = [-150, 40, 0.1; 220, 300, -10];
Np = size(B, 1);
% view angle
viewangle = 0;

% image grid
h = 32.0;
Xgrid = (-8:8).*h;
Ygrid = (-8:8).*h;
Zgrid = (-2:2).*h;
Nx = size(Xgrid, 2);
Ny = size(Ygrid, 2);
Nz = size(Zgrid, 2);
% Lxy is the length of AB on xy plane
Lxy = sqrt((B(:,1)-A(:,1)).^2 + (B(:,2)-A(:,2)).^2);
% d is the distance of AB to ISO
d = (A(:,1).*B(:,2)-B(:,1).*A(:,2))./Lxy;
Lmid = sqrt(A(:,1).^2+A(:,2).^2-d.^2);
% Zctg is the ctg(theta_z) = Lxy/Z_AB
Zctg = Lxy./(B(:,3)-A(:,3));
% Z_A is A(:,3);
Z_A = A(:,3);
% angles
theta = atan2(B(:,2)-A(:,2), B(:,1)-A(:,1));

[dt, Vindex] = linesinimage3D(theta, d, Lxy, Lmid, Z_A, Zctg, Xgrid, Ygrid, Zgrid);

% plot to check
i_chk = 2;
a = zeros((Nx-1)*(Ny-1)*(Nz-1)+1,1);
a(Vindex(i_chk,:)) = dt(i_chk,:);
a = reshape(a(1:end-1), Nx-1, Ny-1, Nz-1)+0.1*max(a(:));
% a(7) = 3;
xrange = [-Nx/2+1, Nx/2-1].*h;
yrange = [-Ny/2+1, Ny/2-1].*h;
figure;
for ii = 1:Nz-1
    subplot(1, Nz-1, ii);
    hold on
    imagesc(xrange, yrange, a(:,:,ii)', [0, max(a(:))]);
    colormap(flipud(colormap('hot')));
    plot([A(i_chk, 1), B(i_chk, 1)], [A(i_chk, 2), B(i_chk, 2)]);
    axis equal
    axis([-500, 500, -500, 500]);
    
end
