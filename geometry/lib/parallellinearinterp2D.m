function [inter_alpha, index_1, index_2, cs_view] = parallellinearinterp2D(Nx, Ny, d_h, viewangle, Nc)

phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);

if phase_flag
    dt = tan(pi/2-viewangle);
    cs_view = csc(viewangle);
    Nx_x = Nx;
    Ny_y = Ny;
else
    dt = tan(viewangle);
    cs_view = -sec(viewangle);
    Nx_x = Ny;
    Ny_y = Nx;
end

Np = size(d_h, 1);
t = repmat((-(Ny_y-1)/2:(Ny_y-1)/2).*dt, Np, 1) - repmat(d_h.*cs_view, 1, Ny_y) + 1/2;
index_1 = floor(t);
inter_alpha = t - index_1;

index_1 = index_1 + Nx_x/2;
index_2 = index_1 + 1;
index_1(index_1<=0 | index_1>Nx_x) = nan;
index_2(index_2<=0 | index_2>Nx_x) = nan;

if phase_flag
    index_1 = index_1 + repmat((0:Ny-1).*Nx, Np, 1);
    index_2 = index_2 + repmat((0:Ny-1).*Nx, Np, 1);
else
    index_1 = (index_1-1).*Nx + repmat(1:Nx, Np, 1);
    index_2 = (index_2-1).*Nx + repmat(1:Nx, Np, 1);
end
index_1(isnan(index_1)) = Nc+1;
index_2(isnan(index_2)) = Nc+1;

return