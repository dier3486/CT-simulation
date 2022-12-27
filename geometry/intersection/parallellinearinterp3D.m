function [interpX, interpY, interpZ, cs_view] = parallellinearinterp3D(imagesize, A, B, theta, center)

if nargin<5
    center = [0 0 0];
end

% rotation of A B
Rtheta = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
A = A*Rtheta + center;
B = B*Rtheta + center;
% phase flag
phase_flag = (theta>pi/4 & theta<=pi*3/4) | (theta>pi*5/4 & theta<=pi*7/4);
if phase_flag
    cs_view = csc(theta);
else
    cs_view = -sec(theta);
end
index_phs = phase_flag+1;

% lines in image 3D 
V = B - A;
Vmax = V(:, index_phs);
AonV = A(:, index_phs)./Vmax;

Dxyz = A - V.*AonV;
Dxyz(:, index_phs) = 0;

Vnr = V./Vmax;
% Nmax = max(Nx, Ny);
% gridmax = -(Nmax-1)/2 : (Nmax-1)/2;
% t_n3 = Vnr(:)*gridmax + Dxyz(:);

gridmax = -(max(imagesize)-1)/2 : (max(imagesize)-1)/2;
interpX = Vnr(:, 1)*gridmax + Dxyz(:, 1) + (imagesize(1)+1)/2;
interpY = Vnr(:, 2)*gridmax + Dxyz(:, 2) + (imagesize(2)+1)/2;
interpZ = Vnr(:, 3)*gridmax + Dxyz(:, 3);

end