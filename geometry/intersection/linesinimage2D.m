function [interpX, interpY, Cs] = linesinimage2D(Nx, Ny, A, B)
% insections of lines in grid-cells image, 2D. LI model
% [interpX, interpY, Cs] = linesinimage2D(Nx, Ny, A, B)
% then to get the projection meaturement by this: (in e.g. intersection.m)
% D = sum(interp2(Cimage, interpX, interpY, 'linear', 0).*Cs, 2) .*L;

Na = size(A,1);
Nb = size(B,1);
Nmax = max(Nx, Ny);

V = B - A;
ABcrs = A(:,1).*B(:,2) - A(:,2).*B(:,1);

[~, phase_index] = max(abs(V(:,1:2)), [], 2);
index_s1 = (phase_index-1).*Nb + (1:Nb)';
index_s1a = (phase_index-1).*Na + (1:Na)';
index_s2 = (2-phase_index).*Nb + (1:Nb)';
s1 = phase_index==1;    % x>=y
% s2 = phase_index==2;    % y>x

Vmax = V(index_s1);
Vnr = V./Vmax;
Dxy = ABcrs./Vmax.*(phase_index.*2-3);

gridmax = -(Nmax-1)/2 : (Nmax-1)/2;
dt = Vnr(index_s2);
t_n2 = dt*gridmax + Dxy;

interpX = repmat(gridmax, Nb, 1);
interpY = interpX;
interpX(~s1, :) = t_n2(~s1, :);
interpY(s1, :) = t_n2(s1, :);

interpX = interpX + (Nx+1)/2;
interpY = interpY + (Ny+1)/2;

Wa = gridmax + 1/2 - reshape(A(index_s1a), Nb, 1);
Wa(Wa>1) = 1;  Wa(Wa<0) = 0;
Wb = gridmax + 1/2 - B(index_s1);
Wb(Wb>1) = 1;  Wb(Wb<0) = 0;
Cs = (Wa - Wb)./Vmax;

return