function D = linesinobject(A, B, crossplane, objecttype, Cimage, flag_filledzero)
% D = linesinobject(A, B, objecttype, Cimage)

if nargin<5
    Cimage = [];
end
if nargin<6
    flag_filledzero = false;
end

% is GPU?
GPUonoff = isa(A, 'gpuArray');
if GPUonoff
    Aclass = classUnderlying(A);
else
    Aclass = class(A);
end
Na = size(A, 1);
[Nb, ~, Nview] = size(B);
% Nview = Nview/3;

if ~isempty(crossplane)
    % planes number
    Npln = size(crossplane, 1);
    % An = A \dot n_plane, Bn = B \dot n_plane, 
    An = reshape(permute(A, [1, 3, 2]), [], 3) * crossplane(:, 1:3)';
    Bn = reshape(permute(B, [1, 3, 2]), [], 3) * crossplane(:, 1:3)';
    An = reshape(An, Na, [], Npln);
    Bn = reshape(Bn, Nb, Nview, Npln);
    % Lcp or Rcp = (d - An)/(Bn - An) which is the intercept of the plane on AB (and on length of AB).
    d = reshape(crossplane(:, 4), [1 1 1 Npln]);
    Lcp = (d - An)./ (Bn - An);
    Rcp = Lcp;
    % whether the AB is cutted on left or right or neither.
    s = An == Bn & An >= d;
    Lcp( An > Bn | s) = -inf;
    Rcp( An < Bn | s) = inf;
    % Rcp( An == Bn & An < d) = -inf;  need not
    % reshape
    Lcp = reshape(Lcp, Nb, Npln, Nview);
    Rcp = reshape(Rcp, Nb, Npln, Nview);
else
    Lcp = [];
    Rcp = [];
end
% the Lcp and Rcp will be catted after the L and R.

switch objecttype
    case 'sphere'
        % a sphere is |r|<=1.
        % GPU buffer (GPU memory is expensive)
        if GPUonoff
            L = zeros(Nb, 2, Nview, Aclass, 'gpuArray');
            R = zeros(Nb, 2, Nview, Aclass, 'gpuArray');
        else
            L = zeros(Nb, 2, Nview, Aclass);
            R = zeros(Nb, 2, Nview, Aclass);
        end
        % Lab^2
        L(:, 2, :) = sum((A - B).^2, 2);
        % |AxB|^2
        R(:, 2, :) = pararea2(A, B);
        % d0^2
        R(:, 1, :) = R(:, 2, :)./L(:, 2, :);
        % d1;
        L(:, 1, :) = sqrt(sum(A.^2, 2) - R(:, 1, :));
        % d2
        R(:, 2, :) = sqrt((1 - R(:, 1, :)).*(R(:, 1, :)<1));
        % Lab
        L(:, 2, :) = sqrt(L(:, 2, :));
        % R1
        R(:, 1, :) = (L(:, 1, :)+R(:, 2, :))./L(:, 2, :);
        % L1
        L(:, 1, :) = (L(:, 1, :)-R(:, 2, :))./L(:, 2, :);
        % R2
        R(:, 2, :) = 1;
        % L2
        L(:, 2, :) = 0;
        
        % cross plane
        L = cat(2, L, Lcp);
        R = cat(2, R, Rcp);
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);        
        D = D.*(D>0);

    case 'cylinder'
        % a cylinder is |z|<=1 & x^2+y^2<=1.
        if GPUonoff
            L = zeros(Nb, 3, Nview, Aclass, 'gpuArray');
            R = zeros(Nb, 3, Nview, Aclass, 'gpuArray');
        else
            L = zeros(Nb, 3, Nview, Aclass);
            R = zeros(Nb, 3, Nview, Aclass);
        end
        % Lxy^2
        L(:, 2, :) = sum((A(:,1:2,:) - B(:,1:2,:)).^2, 2);
        % |AxB|_{xy}^2
        R(:, 2, :) = (A(:,1,:).*B(:,2,:) - A(:,2,:).*B(:,1,:)).^2;
        % d0^2
        R(:, 1, :) = R(:, 2, :)./L(:, 2, :);
        % d1;
        L(:, 1, :) = sqrt(sum(A(:,1:2,:).^2, 2) - R(:, 1, :));
        % d2
        R(:, 2, :) = sqrt((1 - R(:, 1, :)).*(R(:, 1, :)<1));
        % Lab
        L(:, 2, :) = sqrt(L(:, 2, :));
        % R1
        R(:, 1, :) = (L(:, 1, :)+R(:, 2, :))./L(:, 2, :);
        % L1
        L(:, 1, :) = (L(:, 1, :)-R(:, 2, :))./L(:, 2, :);
        % R2 L2
        R(:, 3, :) = cast(B(:, 3, :)<A(:, 3, :), Aclass).*2 - 1;
        R(:, 2, :) = (-R(:, 3, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        L(:, 2, :) = (R(:, 3, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));        
        % R3 L3
        R(:, 3, :) = 1;
        L(:, 3, :) = 0;

        % cross plane
        L = cat(2, L, Lcp);
        R = cat(2, R, Rcp);
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
              
    case 'blade'
        % a blade is 0<=z<=1.
        if GPUonoff
            L = zeros(Nb, 2, Nview, Aclass, 'gpuArray');
            R = zeros(Nb, 2, Nview, Aclass, 'gpuArray');
        else
            L = zeros(Nb, 2, Nview, Aclass);
            R = zeros(Nb, 2, Nview, Aclass);
        end
        % R1 L1
        R(:, 2, :) = cast(B(:, 3, :)<A(:, 3, :), Aclass);
        R(:, 1, :) = (1 - R(:, 2, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        L(:, 1, :) = (R(:, 2, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        % R2 L2
        R(:, 2, :) = 1;
        L(:, 2, :) = 0;

        % cross plane
        L = cat(2, L, Lcp);
        R = cat(2, R, Rcp);
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
        
    case 'cube'
        % a cube is |x|<=1 & |y|<=1 & |z|<=1.
        if GPUonoff
            L = zeros(Nb, 4, Nview, Aclass, 'gpuArray');
            R = zeros(Nb, 4, Nview, Aclass, 'gpuArray');
        else
            L = zeros(Nb, 4, Nview, Aclass);
            R = zeros(Nb, 4, Nview, Aclass);
        end
        % s
        R(:, 1:3, :) = cast(B<A, Aclass).*2 - 1;
        % L123, R123
        L(:, 1:3, :) = (R(:, 1:3, :) - A)./(B - A);
        R(:, 1:3, :) = (-R(:, 1:3, :) - A)./(B - A);
        % L4 R4
        R(:, 4, :) = 1;
        L(:, 4, :) = 0;
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
        
    case 'image2D'
        % 2D image is an image copied on z direction
        % it was a stubid idea to support 'views-lines' and 'views-net',
        % anyhow let's do it
        A = reshape(permute(repmat(A, Nb/Na, 1), [1 3 2]), Nb*Nview, []);
        B = reshape(permute(B, [1 3 2]), Nb*Nview, []);
        % call linesinimage2D
        [Nx, Ny] = size(Cimage);
        if ~flag_filledzero
            [interpX, interpY, Cs] = linesinimage2D(Nx, Ny, A, B);
            D = sum(interp2(Cimage, interpY, interpX, 'linear', 0).*Cs, 2);
        else
            [interpX, interpY, Cs] = linesinimage2D(Nx-2, Ny-2, A, B);
            D = sum(interp2(Cimage, interpY+1, interpX+1, 'linear', 0).*Cs, 2);
        end
    
    case {'image3D', 'images'}
        % 3D image is an array of images on z direction
        % for 'views-lines' and 'views-net' we permute the A
        A = reshape(permute(repmat(A, Nb/Na, 1), [1 3 2]), Nb*Nview, []);
        B = reshape(permute(B, [1 3 2]), Nb*Nview, []);
        [Nx, Ny, Nz] = size(Cimage);
        % call linesinimage3D
       if ~flag_filledzero
            [interpX, interpY, interpZ, Cs] = linesinimage3D(Nx, Ny, Nz, A, B);
            D = sum(interp3(Cimage, interpY, interpX, interpZ, 'linear', 0).*Cs, 2);
       else
            [interpX, interpY, interpZ, Cs] = linesinimage3D(Nx-2, Ny-2, Nz-2, A, B);
            D = sum(interp3(Cimage, interpY+1, interpX+1, interpZ+1, 'linear', 0).*Cs, 2);
        end
        
    otherwise
        D = zeros(Nb, Nview, Aclass);
        return 
end
D = reshape(D, Nb, Nview);

end