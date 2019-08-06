% test script for ideal detector position
% on grid image
% 2D

addpath(genpath('../'));

% focalpos_ideal = [0, -550, 0];
% path AB (as A is the focal spot, B is the detector pixels)
B = [100, -500, 0; 100, -500, 0];
A = [-150, 40, 0; 220, 300, 0];
Np = size(B, 1);
% view angle
viewangle = 0;

% image grid
h = 32.0;
Xgrid = -8:8;
Ygrid = -8:8;
Nx = size(Xgrid,2);
Ny = size(Ygrid,2);
% so the voxel is (Nx-1)*(Ny-1)

% L and d is not changed by viewangle
L = sqrt(sum((B-A).^2, 2));
d = (A(:,1).*B(:,2)-B(:,1).*A(:,2))./L;
% AO = sqrt(L.^2-d.^2);
AO = sqrt(A(:,1).^2+A(:,2).^2-d.^2);
% a normalized L shall be employed (not in this test)
% static angles
theta0 = atan2(B(:,2)-A(:,2), B(:,1)-A(:,1));

% + view angle
theta = mod(theta0 + viewangle, pi*2);
% tan and cotan
tan_theta = tan(theta);
cot_theta = cot(theta);
% delta on path AB
delta_x = h./cos(theta);
delta_y = h./sin(theta);
% the intersection points of the grid with the path AB, along the path
tx = delta_x*Xgrid + repmat(AO-d.*tan_theta, 1, Nx);
ty = delta_y*Ygrid + repmat(AO+d.*cot_theta, 1, Ny);
% sort them, 0 is A, L is B
[t_sort, t_I] = sort([tx, ty, zeros(Np, 1), L], 2);
% 'phase' is the order determined by the direction of AB
phase_x = (theta < pi/2) | (theta >= pi*3/2);
phase_y = theta < pi;
% I know theta in [0, 2*pi)
% t_I to image index
Vx = cumsum(t_I<=Nx, 2);
Vy = cumsum((t_I>Nx) & (t_I<=Nx+Ny), 2);
Vx(Vx==0 | Vx>=Nx) = nan;
Vy(Vy==0 | Vy>=Ny) = nan;
Vx(~phase_x, :) = Nx-Vx(~phase_x, :);
Vy(~phase_y, :) = Ny-Vy(~phase_y, :);
Vindex = Vx.*1 + (Vy-1).*(Nx-1);
naninV = isnan(Vindex);
Vindex(naninV) = (Nx-1)*(Ny-1)+1;

dt = [diff(t_sort,1,2), zeros(Np,1)];
dt = dt.*(cumsum(t_I==Nx+Ny+1, 2)-cumsum(t_I==Nx+Ny+2, 2));

% plot to check
i_chk = 1;
a = zeros((Nx-1)*(Ny-1)+1,1);
a(Vindex(i_chk,:)) = dt(i_chk,:);
a = reshape(a(1:end-1), Nx-1, Ny-1)+0.1*max(a(:));
% a(7) = 3;
xrange = [-Nx/2+1, Nx/2-1].*h;
yrange = [-Ny/2+1, Ny/2-1].*h;
figure;
% plot([A(i_chk, 1), B(i_chk, 1)], [A(i_chk, 2), B(i_chk, 2)]);
hold on
imagesc(xrange, yrange, a', [0, max(a(:))]);
colormap(flipud(colormap('hot')));
plot([A(i_chk, 1), B(i_chk, 1)], [A(i_chk, 2), B(i_chk, 2)]);
axis([-500, 500, -500, 500]);
axis equal



