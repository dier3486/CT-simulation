function D = linesinobjectGPU(A, B, objecttype, Cimage)
% D = linesinobject(A, B, objecttype, Cimage), GPU version

if nargin<4
    Cimage = [];
end

% value class
varclass = classUnderlying(A);
% I know the class of A,B and Cimage are same (in single or double), and they are all gpuArray

N = size(A,1);
switch objecttype
    case 'sphere'
        % a sphere is |r|<=1.
        Lab = sqrt(sum((A - B).^2, 2));
        d0 = pararea(A, B)./Lab;
        d1 = sqrt(sum(A.^2, 2) - d0.^2);
        d2 = 1 - d0.^2;
        d2(d2<0) = 0;
        d2 = sqrt(d2);

        L = [(d1 - d2)./Lab, zeros(N, 1)];
        R = [(d1 + d2)./Lab, ones(N, 1)];
        D = multiinsect(L, R);

    case 'cylinder'
        % a cylinder is |z|<=1 & x^2+y^2<=1.
        Lxy = sqrt(sum((A(:,1:2) - B(:,1:2)).^2, 2));
        S = abs(A(:,1).*B(:,2) - A(:,2).*B(:,1));
        d0 = S./Lxy;
        d1 = sqrt(sum(A(:,1:2).^2, 2) - d0.^2);
        d2 = 1 - d0.^2;
        d2(d2<0) = 0;
        d2 = sqrt(d2);
        L1 = (d1-d2)./Lxy;
        R1 = (d1+d2)./Lxy;

        LR2 = [-1 - A(:,3), 1 - A(:,3)]./repmat(B(:,3) - A(:,3), 1, 2);
        sn = B(:,3)<A(:,3);
        LR2(sn, :) = fliplr(LR2(sn, :)); 

        L = [L1, LR2(:,1), zeros(N, 1)];
        R = [R1, LR2(:,2), ones(N, 1)];
        D = multiinsect(L, R);
        
    case 'blade'
        % a blade is 0<=z<=1.
        LR = [-A(:,3), 1 - A(:,3)]./repmat(B(:,3) - A(:,3), 1, 2);
        sn = B(:,3)<A(:,3);
        LR(sn, :) = fliplr(LR(sn, :)); 

        L = [LR(:,1), zeros(N, 1)];
        R = [LR(:,2), ones(N, 1)];
        D = multiinsect(L, R);
        
    case 'cube'
        % a cube is |x|<=1 & |y|<=1 & |z|<=1.
        L1 = (-1 - A)./(B - A);
        R1 = (1 - A)./(B - A);
        sn = B<A;
        L(~sn) = L1(~sn);
        L(sn) = R1(sn);
        R(~sn) = R1(~sn);
        R(sn) = L1(sn);
        L = [L, zeros(N, 1)];
        R = [R, ones(N,1)];
        D = multiinsect(L, R);
        
    case 'image2D'
        % 2D image is an image copied on z direction
        % We strongly suggest to call projectioninimage.m in projection
        % simulations but not this for performance
        % grid
        [Nx, Ny] = size(Cimage);
        Xgrid = -Nx/2:Nx/2;
        Ygrid = -Ny/2:Ny/2;        
        % Lxy is the length of AB on xy plane
        Lxy = sqrt((B(:,1)-A(:,1)).^2 + (B(:,2)-A(:,2)).^2);
        % d is the distance from ISO to AB
        d = (A(:,2).*B(:,1)-A(:,1).*B(:,2))./Lxy;
        Lmid = sqrt(A(:,1).^2+A(:,2).^2-d.^2);
        % angles
        theta = atan2(B(:,2)-A(:,2), B(:,1)-A(:,1));
        % call 2D projection function
        [dt, Vindex] = linesinimage2D(theta, d, Lxy, Lmid, Xgrid, Ygrid);
        Cimage = [Cimage(:); 0];
        D = sum(dt.*Cimage(Vindex), 2)./Lxy;
        % here we used lines' insection method as a projection
    
    case {'image3D', 'images'}
        % 3D image is an array of images on z direction
        % We strongly suggest to call projectioninimage.m in projection
        % simulations but not this for performance
        % grid
        [Nx, Ny, Nz] = size(Cimage);
        Xgrid = -Nx/2:Nx/2;
        Ygrid = -Ny/2:Ny/2;
        Zgrid = -Nz/2:Nz/2;
        % Lxy is the length of AB on xy plane
        Lxy = sqrt((B(:,1)-A(:,1)).^2 + (B(:,2)-A(:,2)).^2);
%         % d is the distance of AB to ISO
%         d = (A(:,1).*B(:,2)-B(:,1).*A(:,2))./Lxy;
        % d is the distance from ISO to AB
        d = (A(:,2).*B(:,1)-A(:,1).*B(:,2))./Lxy;
        Lmid = sqrt(A(:,1).^2+A(:,2).^2-d.^2);
        % Zctg is the ctg(theta_z) = Lxy/Z_AB
        Zctg = Lxy./(B(:,3)-A(:,3));
        % Z_A is A(:,3);
        Z_A = A(:,3);
        % angles
        theta = atan2(B(:,2)-A(:,2), B(:,1)-A(:,1));
        % call 3D projection function
        [dt, Vindex] = linesinimage3D(theta, d, Lxy, Lmid, Z_A, Zctg, Xgrid, Ygrid, Zgrid);
        Cimage = [Cimage(:); 0];
        D = sum(dt.*Cimage(Vindex), 2)./Lxy;
        % here we used lines' insection method as a projection
        
    otherwise
        D = zeros(N, 1);
        return 
end


return