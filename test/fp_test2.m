Cimage = phantom().*100;
Oimg = [0 0 0];

[Nx, Ny] = size(Cimage);

Np = 500;
A = [0 -200 0];
B = [linspace(-200, 200, Np)' ones(Np,1).*200 zeros(Np,1)];
h = 1;

Nview = 1152;
views = linspace(0, pi*2*(Nview-1)/Nview, Nview);

maxview = 32;


% go
Cimage_ext = [Cimage(:); 0];
Next = Nx*Ny+1;
% Lxy is the length of AB on xy plane
Lxy = sqrt((B(:,1)-A(:,1)).^2 + (B(:,2)-A(:,2)).^2);
% Lo is the length from ISO to A (on xy plane)
Lo = sqrt(A(:,1).^2 + A(:,2).^2);

% d is the distance from ISO to AB
d = (A(:,2).*B(:,1)-A(:,1).*B(:,2))./Lxy;
Lmid = sqrt(A(:,1).^2+A(:,2).^2-d.^2);
        
% Xgrid = -(Nx-1)/2:(Nx-1)/2;
% Ygrid = -(Ny-1)/2:(Ny-1)/2; 

Rt = 141;
t = 1.0;
tgrid = [fliplr(0:-t:-Rt) (t:t:Rt)]' + Lo;
Nt = length(tgrid);

Psample = reshape(reshape((B - A)./Lxy, [], 1)*tgrid', Np, 3, Nt) + A;
Psample = reshape(permute(Psample, [1 3 2]), [], 3);

Rp2 = sum(Psample.^2, 2);
Sp = Rp2<Rt^2;
P0 = Psample(Sp, :);
Nsp = sum(Sp);

Nv = maxview;
v1 = 1;
v2 = v1+Nv-1;
viewangle = views(v1:v2);
tiltangle = zeros(Nv, 1);

tic;
viewangle = reshape(viewangle, 1, 1, []);
tiltangle = reshape(tiltangle, 1, 1, []);
Mv = [cos(viewangle)    sin(viewangle).*cos(tiltangle)   sin(viewangle).*sin(tiltangle);
    -sin(viewangle)    cos(viewangle).*cos(tiltangle)   cos(viewangle).*sin(tiltangle);
     zeros(1,1,Nv)     -sin(tiltangle)                   cos(tiltangle)                ];
Mv = reshape(Mv, 3, []);

Pv = P0*Mv - 0.5;
Pv = reshape(permute(reshape(Pv, Nsp, 3, Nv), [1 3 2]), [], 3);

index_xy = floor(Pv(:, [1 2]));
alpha_xy = Pv(:, [1 2]) - index_xy;
index_xy = index_xy + [Nx/2+1 Ny/2];
sout = any((index_xy<[1 0]) | (index_xy>[Nx-1 Ny-2]), 2);

index_xy1 = index_xy*[1; Nx];
index_intp = [index_xy1  index_xy1+1  index_xy1+Nx  index_xy1+Nx+1];
index_intp(sout, :) = Next;
alpha_intp = [(1-alpha_xy(:,1)).*(1-alpha_xy(:,2))  alpha_xy(:,1).*(1-alpha_xy(:,2))  ...
    (1-alpha_xy(:,1)).*alpha_xy(:,2)  alpha_xy(:,1).*alpha_xy(:,2)];

Dintp = sum(Cimage_ext(index_intp).*alpha_intp, 2);
Dintp = reshape(Dintp, Nsp, Nv);

Dsample = zeros(Np*Nt, Nv);
Dsample(Sp, :) = Dintp;
Dsample = reshape(Dsample, Np, Nt, Nv);

Dv = squeeze(sum(Dsample, 2).*t);
toc;

tic;
% 0 & 1
D0 = zeros(Np, Nv);
D1 = zeros(Np, Nv);
d_h = d./h;
Xgrid = -Nx/2:Nx/2;
Ygrid = -Ny/2:Ny/2;
for iview = v1:v2
    Mrot = [cos(views(iview))  sin(views(iview))    0;
           -sin(views(iview))  cos(views(iview))    0;
            0                  0                    1];
    % angles
    Arot = A*Mrot;
    Brot = B*Mrot;
    theta = atan2(Brot(:,2)-Arot(:,2), Brot(:,1)-Arot(:,1));
    % 0
    [dt, Vindex] = linesinimage2D(theta, d, Lxy, Lmid, Xgrid, Ygrid);
    D0(:, iview) = sum(dt.*Cimage_ext(Vindex), 2);
    % 1
    [inter_alpha, index_1, index_2, cs_vangle] = linesinimageLI2D(Nx, Ny, d_h, theta);
    D1(:, iview) = sum(Cimage_ext(index_1).*(1-inter_alpha) + Cimage_ext(index_2).*inter_alpha, 2).*cs_vangle;
end
toc;
