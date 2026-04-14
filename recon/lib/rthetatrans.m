function raw = rthetatrans(img, varargin)
% r-theta transform of the img
% raw = rthetatrans(img, center, Ntheta, d, flag_even);
% or, raw = rthetatrans(img);

% default inputs
%               center  Ntheta   d        R,       flag_even
defaultinput = {[0 0],  192,     1.0,     inf,     true      };
% input coeffients
[center, Ntheta, d, R, flag_even] = cleaninputarg(defaultinput, varargin{:});

% image size
[Ny, Nx, Nrow] = size(img);

% rep center
if size(center, 1) == 1
    center = repmat(center, Nrow, 1);
end

% XY grid of image
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(:, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(:, 2);
% theta grid
Vb = (0:Ntheta-1).*pi/Ntheta - pi/2;
% R grid
R = min(R, max(Nx, Ny)*sqrt(2)/2);
if flag_even
    Nb = ceil(R/d)*2;
else
	Nb = ceil(R/d)*2 + 1;
end
Rb = (-(Nb-1)/2 : (Nb-1)/2).*d;

% interp target 
Xb = Rb(:) * cos(Vb);
Yb = Rb(:) * sin(Vb);

% interp 2D
raw = zeros(Nb, Ntheta, Nrow, 'like', img);
for irow = 1:Nrow
	raw(:, :, irow) = interp2(Xa(irow, :), Ya(irow, :), img(:, :, irow), Xb, Yb, 'linear', 0);
end
% done

end