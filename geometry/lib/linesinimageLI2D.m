function [inter_alpha, index_1, index_2, cs_vangle] = linesinimageLI2D(Nx, Ny, d_h, theta, Next)
% insections of lines in grid-cells image, 2D. linear interp
% [dt, Vindex] = linesinimage2D(theta, d, L, AO, Xgrid, Ygrid);
% then D = sum(dt.*Cimage(Vindex), 2); 
% remember to add a 0 after Cimage that Cimage = [Cimage(:); 0]
% for infinite lines, use:
% [dt, Vindex] = linesinimage2D(theta, d, [], 0, Xgrid, Ygrid);

if nargin<5
    Next = Nx*Ny+1;
end

% mod theta with 2pi
theta = mod(theta(:), pi*2);
% I know theta = theta(:); d = d(:);

% ini
Np = size(theta, 1);
N = max(Nx, Ny);
cs_vangle = zeros(Np, 1);
inter_alpha = zeros(Np, N);
index_1 = ones(Np, N).*Next;
index_2 = ones(Np, N).*Next;

% loop 2 phase
phase_flag = (theta>pi/4 & theta<=pi*3/4) | (theta>pi*5/4 & theta<=pi*7/4);
% phase 1

for iphase = [1 0]
    if iphase
        vangle = theta(phase_flag);
        d_i = d_h(phase_flag);
        dt_i = tan(pi/2-vangle);
        cs_vangle_i = csc(vangle);
        cs_vangle(phase_flag) = cs_vangle_i;
        Nx_x = Nx;
        Ny_y = Ny;
    else
        vangle = theta(~phase_flag);
        d_i = d_h(~phase_flag);
        dt_i = tan(vangle);
        cs_vangle_i = -sec(vangle);
        cs_vangle(~phase_flag) = cs_vangle_i;
        Nx_x = Ny;
        Ny_y = Nx;
    end
    
    t = dt_i*(-(Ny_y-1)/2:(Ny_y-1)/2) - repmat(d_i.*cs_vangle_i, 1, Ny_y) + 1/2;
    index_1_i = floor(t);
    inter_alpha_i = t - index_1_i;

    index_1_i = index_1_i + Nx_x/2;
    index_2_i = index_1_i + 1;
    s1 = index_1_i<=0 | index_1_i>Nx_x;
    s2 = index_2_i<=0 | index_2_i>Nx_x;
    if iphase
        index_1_i = index_1_i + (0:Ny-1).*Nx;
        index_2_i = index_2_i + (0:Ny-1).*Nx;
    else
        index_1_i = (index_1_i-1).*Nx + (1:Nx);
        index_2_i = (index_2_i-1).*Nx + (1:Nx);
    end
    index_1_i(s1) = Next;
    index_2_i(s2) = Next;
    
    if iphase
        inter_alpha(phase_flag, 1:Ny) = inter_alpha_i;
        index_1(phase_flag, 1:Ny) = index_1_i;
        index_2(phase_flag, 1:Ny) = index_2_i;
    else
        inter_alpha(~phase_flag, 1:Nx) = inter_alpha_i;
        index_1(~phase_flag, 1:Nx) = index_1_i;
        index_2(~phase_flag, 1:Nx) = index_2_i;
    end
    
end

return