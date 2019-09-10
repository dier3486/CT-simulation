function [dt, Vindex] = linesinimage3D(theta, d, Lxy, Lmid, Z_A, Zctg, Xgrid, Ygrid, Zgrid)
% insections of lines in grid-cells image, 3D.
% [dt, Vindex] = linesinimage2D(theta, d, Lxy, Lmid, Z_A, Zctg, Xgrid, Ygrid, Zgrid);
% then D = sum(dt.*Cimage(Vindex), 2).*L./Lxy; 
% remember to add a 0 after Cimage that Cimage = [Cimage(:); 0]

% the numbers
N = size(theta, 1);
Nx = size(Xgrid(:), 1);
Ny = size(Ygrid(:), 1);
Nz = size(Zgrid(:), 1);
% mod theta with 2pi
theta = mod(theta, pi*2);
% % tan and cotan
% tan_theta = tan(theta);
% cot_theta = cot(theta);
% delta on path AB
h = 1.0;
% delta_x = h.*sec(theta);    % h./cos(theta)
% delta_y = h.*csc(theta);    % h./sin(theta)
delta_z = h.*Zctg;
% the intersection points of the grid with the path AB, along the path
% tx = delta_x*Xgrid + repmat(Lmid + d.*tan_theta, 1, Nx);
% ty = delta_y*Ygrid + repmat(Lmid - d.*cot_theta, 1, Ny);
tx = (Xgrid.*h + d.*sin(theta)).*sec(theta) + Lmid;
ty = (Ygrid.*h - d.*cos(theta)).*csc(theta) + Lmid;
if size(Z_A, 1)>1
    tz = (repmat(Zgrid, N, 1)-repmat(Z_A, 1, Nz)).*repmat(delta_z, 1,Nz);
else
    tz = (repmat(Zgrid, N, 1)-Z_A).*repmat(delta_z, 1,Nz);
end
% sort them, 0 is A, L is B
if isempty(Lxy)
    [t_sort, t_I] = sort([tx, ty, tz], 2);
else
    [t_sort, t_I] = sort([tx, ty, tz, zeros(N, 1), Lxy], 2);
end
% 'phase' is the order determined by the direction of AB
phase_x = (theta < pi/2) | (theta >= pi*3/2);
phase_y = theta < pi;
phase_z = Zctg > 0;
% I know theta in [0, 2*pi)
% t_I to image index
Vx = cumsum(t_I<=Nx, 2);
Vy = cumsum((t_I>Nx) & (t_I<=Nx+Ny), 2);
Vz = cumsum((t_I>Nx+Ny) & (t_I<=Nx+Ny+Nz), 2);
Vx(Vx==0 | Vx>=Nx) = nan;
Vy(Vy==0 | Vy>=Ny) = nan;
Vz(Vz==0 | Vz>=Nz) = nan;
Vx(~phase_x, :) = Nx-Vx(~phase_x, :);
Vy(~phase_y, :) = Ny-Vy(~phase_y, :);
Vz(~phase_z, :) = Nz-Vz(~phase_z, :);
% Vindex
Vindex = Vx.*1 + (Vy-1).*(Nx-1) + (Vz-1).*((Nx-1)*(Ny-1));
naninV = isnan(Vindex);
Vindex(naninV) = (Nx-1)*(Ny-1)*(Nz-1)+1;
% dt
dt = [diff(t_sort,1,2), zeros(N,1)];
dt(~isfinite(dt)) = 0;
if ~isempty(Lxy)
    dt = dt.*(cumsum(t_I==Nx+Ny+Nz+1, 2) - cumsum(t_I==Nx+Ny+Nz+2, 2));
end
% D = sum(dt.*Cimage(Vindex), 2);

return