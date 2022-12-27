function [interpX, interpY, interpZ, Cs] = linesinimage3D(Nx, Ny, Nz, A, B)
% insections of lines in grid-cells image, 3D. LI model
% [interpX, interpY, interpZ, Cs] = linesinimage3D(Nx, Ny, Nz, A, B);
% then to get the projection meaturement by this: (in e.g. intersection.m)
% D = sum(interp3(Cimage, interpX, interpY, interpZ, 'linear', 0).*Cs, 2) .*L;

Na = size(A,1);
Nb = size(B,1);
% I know Na=1 or Na=Nb
Nmax = max([Nx, Ny, Nz]);
% the returned interpXYZ will in size (Nb, Nmax)

V = B - A;
[~, phase_index] = max(abs(V), [], 2);
index_s1 = (phase_index-1).*Nb + (1:Nb)';
index_s1a = (phase_index-1).*Na + (1:Na)';

Vmax = V(index_s1);
AonV = reshape(A(index_s1a), Nb, 1)./Vmax;
Dxyz = A - V.*AonV;

Vnr = V./Vmax;
gridmax = -(Nmax-1)/2 : (Nmax-1)/2;
t_n3 = Vnr(:)*gridmax + Dxyz(:);

interpX = t_n3(1:Nb, :) + (Nx+1)/2;
interpY = t_n3(Nb+1:Nb*2, :) + (Ny+1)/2;
interpZ = t_n3(Nb*2+1:Nb*3, :) + (Nz+1)/2;

Wa = gridmax + 1/2 - reshape(A(index_s1a), Nb, 1);
Wa(Wa>1) = 1;  Wa(Wa<0) = 0;
Wb = gridmax + 1/2 - B(index_s1);
Wb(Wb>1) = 1;  Wb(Wb<0) = 0;
% W = Wa - Wb;
% Cs = 1./Vmax;
Cs = (Wa - Wb)./Vmax;

return