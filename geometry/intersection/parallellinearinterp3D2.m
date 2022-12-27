function [interpX, interpY, interpZ, cs_view] = parallellinearinterp3D2(imagesize, A, V, viewangle, Zgrid, centerXY_h)

if nargin<6
    centerXY_h = [0 0];
end
if (viewangle>pi*2)
    error('Please mod the viewangle by 2*pi!');
end

Nz = length(Zgrid);
Np = size(A, 1);
phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);
if phase_flag
    cs_view = csc(viewangle);
else
    cs_view = -sec(viewangle);
end

Rtheta = [cos(viewangle) sin(viewangle); -sin(viewangle) cos(viewangle)];
A = A*Rtheta + centerXY_h;
V = V*Rtheta;

index_p = phase_flag+1;
index_n = ~phase_flag+1;
invVmax = 1./V(:, index_p);
V(:, index_p) = 1;
V(:, index_n) = V(:, index_n).*invVmax;

Vnr = zeros(Np*Nz, 3, 'like', A);
Vnr(:, [1 2]) = repmat(V, Nz, 1);
Vnr(:, 3) = reshape(invVmax*Zgrid, [], 1);

Dxyz = zeros(Np*Nz, 3, 'like', A);
% Dxyz(:, index_p) = 0;
Dxyz(:, index_n) = repmat(A(:, index_n)-V(:, index_n).*A(:, index_p), Nz, 1);
Dxyz(:, 3) = -reshape((invVmax.*A(:, index_p))*Zgrid, [], 1);

gridmax = -(max(imagesize)-1)/2 : (max(imagesize)-1)/2;
interpX = Vnr(:, 1)*gridmax + Dxyz(:, 1) + (imagesize(1)+1)/2;
interpY = Vnr(:, 2)*gridmax + Dxyz(:, 2) + (imagesize(2)+1)/2;
interpZ = Vnr(:, 3)*gridmax + Dxyz(:, 3);

end