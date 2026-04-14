function u = TVfunG(f, u, b, d, ~, Cl, balance, Ndim, TVopts, BoundCond)

if nargin < 10
    BoundCond = TVboundarycond(TVopts);
end

sizeU = size(u);

u1 = 0;
for idim = 1 : Ndim
    if sizeU(idim) == 1
        continue;
    end
    % dimorder
    dimorder = circshift(1:3, -idim+1);
    % re-order the u
    u = permute(u, dimorder);
    u_idim = u([2:end end], :, :) + u([1 1:end-1], :, :);
    % boundary condition left
    switch BoundCond(idim, 1)
        case 0
            % skip
        case 1 % clamp
            % did
        case 2 % linear extrap
            u_idim(1, :, :) = u(1, :, :).*2;
        case 3 % mirror
            u_idim(1, :, :) = u(2, :, :).*2;
        otherwise
            error(0);
    end
    % boundary condition right
    switch BoundCond(idim, 2)
        case 0
            % skip
        case 1 % clamp
            % did
        case 2 % linear extrap
            u_idim(end, :, :) = u(end, :, :).*2;
        case 3 % mirror
            u_idim(end, :, :) = u(end-1, :, :).*2;
        otherwise
            error(0);
    end

    % add to u1
    u1 = u1 + ipermute(u_idim, dimorder).*balance(idim);
    % inverse the order
    u = ipermute(u, dimorder);

    % re-order the b,d to calculate d(b,d)
    b = permute(b, [dimorder, 4]);
    d = permute(d, [dimorder, 4]);
    b_idim = b(1:end,:,:, idim) - b([1 1:end-1],:,:, idim);
    d_idim = d(1:end,:,:, idim) - d([1 1:end-1],:,:, idim);
    % boundary condition left
    switch BoundCond(idim, 1)
        case 0
            % skip
        case 1 % clamp
            b_idim(1, :, :) = b(1,:,:, idim);
            d_idim(1, :, :) = d(1,:,:, idim);
        case 2 % linear extrap
            % did
        case 3 % mirror
            b_idim(1, :, :) = b(1,:,:, idim).*2;
            d_idim(1, :, :) = d(1,:,:, idim).*2;
        otherwise
            error(0);
    end
    u1 = u1 + ipermute(b_idim - d_idim, dimorder).*balance(idim);
    % inverse the order
    b = ipermute(b, [dimorder, 4]);
    d = ipermute(d, [dimorder, 4]);
end

% fix boundary of u
leftFix = BoundCond(:, 1) == 0;
rightFix = BoundCond(:, 2) == 0;
% new u
u(1+leftFix(1) : end-rightFix(1), 1+leftFix(2) : end-rightFix(2), 1+leftFix(3) : end-rightFix(3)) = ...
    u1(1+leftFix(1) : end-rightFix(1), 1+leftFix(2) : end-rightFix(2), 1+leftFix(3) : end-rightFix(3)).*(Cl/(1+Cl*Ndim*2)) + ...
    f(1+leftFix(1) : end-rightFix(1), 1+leftFix(2) : end-rightFix(2), 1+leftFix(3) : end-rightFix(3)).*(1/(1+Cl*Ndim*2));

end