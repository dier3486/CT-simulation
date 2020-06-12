function [A, L, indexrange] = pinleftoptmatrix(cs, Npixel, Nslice, Nview, Npixelpermod, alpha)
% subroutine used in geometric cali

A = cell(Nslice, 1);
indexrange = zeros(Nslice, 2);
for islice = 1:Nslice
    pindex = [cs(islice, :).pindex];
    vindex = [cs(islice, :).vindex];
    A{islice} = sparse(vindex, pindex, double([cs(islice, :).dp]), Nview, Npixel);
    % cut edge
    edge1 = (ceil(min(pindex)/Npixelpermod) + 2)*Npixelpermod + 1;
    edge2 = (floor(max(pindex)/Npixelpermod) - 2)*Npixelpermod;
    indexrange(islice, :) = [edge1, edge2];
%     A{islice} = A{islice}(:, p1:p2);
end

% dgA = zeros(Npixel, Nslice);
% L = spdiags(ones(Npixel, 1)*[-1 2 1], [-1 0 1], Npixel, Npixel);
Nmod = Npixel/Npixelpermod;
Lmod = cell(1, Nmod);
Lvalue = repmat([-1/2 1 -1/2], Npixelpermod, 1);
Lvalue(1, 2) = 1/2;  Lvalue(end, 2) = 1/2;  
Lmod(:) = {spdiags(Lvalue, [-1 0 1], Npixelpermod, Npixelpermod)};
Lall = spdiags(repmat([-1/2 1 -1/2], Npixel, 1), [-1 0 1], Npixel, Npixel);
Lall(1, 1) = 1/2;  Lall(end, end) = 1/2;
alpha_mod = alpha(1)-alpha(2);  alpha_all = alpha(2);
L = blkdiag(Lmod{:}).*alpha_mod + Lall.*alpha_all;
end