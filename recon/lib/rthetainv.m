function img = rthetainv(raw, varargin)
% inverse r-theta transform
% img = rthetainv(raw, center, imgsize, d, flag_fill);
% or, img = rthetainv(raw);

% raw size
[Nb, Ntheta, Nrow] = size(raw);

% default inputs
%               center           imgsize   d        flag_fill
defaultinput = {zeros(Nrow, 2),  512,      1.0,     true      };
% input coeffients
[center, imgsize, d, flag_fill] = cleaninputarg(defaultinput, varargin{:});
% rep center
if size(center, 1) == 1
    center = repmat(center, Nrow, 1);
end

% XY grid of image(s)
Nx = imgsize(1);
Ny = imgsize(end);
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(:, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(:, 2);
% theta grid
if flag_fill
    Vb = (0:Ntheta).*pi/Ntheta - pi/2;
else
    Vb = (0:Ntheta-1).*pi/(Ntheta-1) - pi/2;
end
% R grid
Rb = (-(Nb-1)/2 : (Nb-1)/2).*d;

% ini img
img = zeros(Nx, Ny, Nrow, 'like', raw);
% loop the rows
for irow = 1:Nrow
    % interp target
    [ndXa, ndYa] = ndgrid(Xa(irow, :), Ya(irow, :));
    Va = fillmissing(atan(ndYa./ndXa), 'constant', 0);
    Ra = sqrt(ndYa.^2 + ndXa.^2);
    Ra(ndXa<0) = -Ra(ndXa<0);
    if flag_fill
        img(:, :, irow) = interp2(Vb, Rb, [raw(:, :, irow) flipud(raw(:, 1, irow))], Va, Ra);
    else
        img(:, :, irow) = interp2(Vb, Rb, raw(:, :, irow), Va, Ra);
    end
    
end


% done

end