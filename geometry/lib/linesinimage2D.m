function [dt, Vindex] = linesinimage2D(theta, d, L, AO, Xgrid, Ygrid)
% insections of lines in grid-cells image, 2D. mosaic modle
% [dt, Vindex] = linesinimage2D(theta, d, L, AO, Xgrid, Ygrid);
% then D = sum(dt.*Cimage(Vindex), 2); 
% remember to add a 0 after Cimage that Cimage = [Cimage(:); 0]
% for infinite lines, use:
% [dt, Vindex] = linesinimage2D(theta, d, [], 0, Xgrid, Ygrid);

% the numbers
N = size(theta, 1);
Nx = size(Xgrid(:),1);
Ny = size(Ygrid(:),1);
% set h = 1
h = 1.0;
% mod theta with 2pi
theta = mod(theta, pi*2);
% % delta on path AB
% delta_x = h.*sec(theta);    % h./cos(theta)
% delta_y = h.*csc(theta);    % h./sin(theta)

% the intersection points of the grid with the path AB, along the path
% tx = delta_x*Xgrid + repmat(AO + d.*tan_theta, 1, Nx);
% ty = delta_y*Ygrid + repmat(AO - d.*cot_theta, 1, Ny);
tx = (Xgrid.*h + d.*sin(theta)).*sec(theta) + AO;
ty = (Ygrid.*h - d.*cos(theta)).*csc(theta) + AO;
% sort them, 0 is A, L is B
if isempty(L)
    [t_sort, t_I] = sort([tx, ty], 2);
else
    [t_sort, t_I] = sort([tx, ty, zeros(N, 1), L], 2);
%     t_sort = [tx, ty, zeros(N, 1), L];
%     t_I = repmat(1:Nx+Ny+2, N, 1);
end
% 'phase' is the order determined by the direction of AB
phase_x = (theta <= pi/2) | (theta > pi*3/2);
phase_y = theta <= pi;
% I know theta in [0, 2*pi)
% t_I to image index
Vx = cumsum(t_I<=Nx, 2);
Vy = cumsum((t_I>Nx) & (t_I<=Nx+Ny), 2);
Vx(Vx==0 | Vx>=Nx) = nan;
Vy(Vy==0 | Vy>=Ny) = nan;
Vx(~phase_x, :) = Nx-Vx(~phase_x, :);
Vy(~phase_y, :) = Ny-Vy(~phase_y, :);
% Vindex
Vindex = Vx.*1 + (Vy-1).*(Nx-1);
naninV = isnan(Vindex);
Vindex(naninV) = (Nx-1)*(Ny-1)+1;
% dt
dt = [diff(t_sort,1,2), zeros(N,1)];
dt(~isfinite(dt)) = 0;
if ~isempty(L)
    dt = dt.*(cumsum(t_I==Nx+Ny+1, 2) - cumsum(t_I==Nx+Ny+2, 2));
end
% D = sum(dt.*Cimage(Vindex), 2);

return