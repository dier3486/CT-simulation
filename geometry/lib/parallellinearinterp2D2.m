function [interpX, interpY, cs_view] = parallellinearinterp2D2(Nx, Ny, d_h, viewangle, centerXY_h)

if nargin<5
    centerXY_h = [0 0];
end

phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);
% direct_flag = ~(viewangle>pi/4 & viewangle<=pi*5/4);

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
d_h = d_h + centerXY_h*[-sin(viewangle); cos(viewangle)];
t = repmat((-(Ny_y-1)/2:(Ny_y-1)/2).*dt, Np, 1) - repmat(d_h.*cs_view, 1, Ny_y) + 1/2;
if phase_flag
    interpX = t + Nx_x/2;
    interpY = repmat(1:Ny_y, Np, 1);
else
    interpY = t + Ny_y/2;
    interpX = repmat(1:Nx_x, Np, 1);
end

end