function [d1, b1] = TVfundbyub(u, b0, lambda, Ndim, BoundCond)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)
% d_y^{k+1} = shrink(D_yu^{k+1}+b_y^{k}, 1/\lambda)
% d_z^{k+1} = shrink(D_zu^{k+1}+b_z^{k}, 1/\lambda)

if nargin < 5
    BoundCond = ones(3,2);
end

sizeU = size(u);
du = b0.*0;
dimorder = 1:3;
for idim = 1:Ndim
    if sizeU(idim) == 1
        continue;
    end
    % re-order the u to calculate du
    dimorder(1:Ndim) = circshift(1:Ndim, -idim+1);
    u = permute(u, dimorder);
    du_idim = u([2:end end], :, :) - u([1:end-1  end-1], :,:);
    % boundary condition (right)
    switch BoundCond(idim, 2)
        case {0, 1} % clamp
            du_idim(end,:,:) = 0;
        case 2 % linear extrap
            % did
        case 3 % mirror
            du_idim(end,:,:) = -du_idim(end,:,:);
        case 4 % odd extrap
            du_idim(end,:,:) = (u(end, :, :) - u(end-2, :, :))./2;
        otherwise
            error(0);
    end
        
    % permute back
    du(:,:,:, idim) = ipermute(du_idim, dimorder);
    u = ipermute(u, dimorder);
end

b0 = b0 + du;
absb0 = abs(b0);
absb0 = absb0 + (absb0<eps);   % to replace fillmissing
d1 = max(absb0-1./lambda, 0).*(b0./absb0);

b1 = b0 - d1;

end