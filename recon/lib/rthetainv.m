function img = rthetainv(raw, varargin)
% inverse r-theta transform
% img = rthetainv(raw, center, imgsize, d, flag_fill);
% or, img = rthetainv(raw);

% raw size
[Nb, Ntheta, Nrow] = size(raw);

% default inputs
%               center                        imgsize   d        flag_fill
defaultinput = {zeros(Nrow, 2, 'like', raw),  512,      1.0,     true      };
% input coeffients
[center, imgsize, d, flag_fill] = cleaninputarg(defaultinput, varargin{:});
% rep center
if size(center, 1) == 1
    center = repmat(center, Nrow, 1);
end

% XY grid of image(s)
Nx = imgsize(1);
Ny = imgsize(end); % I know the 'end' is 1 or 2
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(:, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(:, 2);
% theta grid
if flag_fill
    Vb = (0:Ntheta).*pi/Ntheta - pi/2;
else
    Vb = (0:Ntheta-1).*pi/(Ntheta-1) - pi/2;
end
Vb = cast(Vb, 'like', raw);
% R grid
Rb = cast((-(Nb-1)/2 : (Nb-1)/2).*d, 'like', raw);

% ini img
img = zeros(Ny, Nx, Nrow, 'like', raw);
% loop the rows
for irow = 1:Nrow
    % interp target
    [msXa, msYa] = meshgrid(Xa(irow, :), Ya(irow, :));
    Va = fillmissing(atan(msYa./msXa), 'constant', 0);
    Ra = sqrt(msYa.^2 + msXa.^2);
    Ra(msXa<0) = -Ra(msXa<0);
    if flag_fill
        img(:, :, irow) = interp2(Vb, Rb, [raw(:, :, irow) flipud(raw(:, 1, irow))], Va, Ra);
    else
        img(:, :, irow) = interp2(Vb, Rb, raw(:, :, irow), Va, Ra);
    end
    
end

% done
end