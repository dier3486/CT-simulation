function raw = rthetatrans(img, center, Ntheta, d, flag_even)
% r-theta transform of the img
% raw = rthetatrans(img, center, Ntheta, d, flag_even);
% or, raw = rthetatrans(img); for lazy guy.

% image size
[Nx, Ny, Nrow] = size(img);

% default prm
if nargin<2 || isempty(center)
	center = zeros(Nrow, 2);
end
if nargin<3
	Ntheta = 180;
end
if nargin<4
	d = 1.0;
end
if nargin<5
	flag_even = true;
end

% XY grid of image
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(:, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(:, 2);
% theta grid
Vb = (0:Ntheta-1).*pi/Ntheta - pi/2;
% R grid
Nb = ceil(max(Nx, Ny)/d*sqrt(2)/2)*2;
if ~flag_even
	Nb = Nb+1;
end
Rb = (-(Nb-1)/2 : (Nb-1)/2).*d;

% interp target 
Xb = Rb(:) * cos(Vb);
Yb = Rb(:) * sin(Vb);

% interp 2D
raw = zeros(Nb, Ntheta, Nrow, 'like', img);
for irow = 1:Nrow
	raw(:, :, irow) = interp2(Xa(irow, :), Ya(irow, :), img(:, :, irow).', Xb, Yb);
end
% done

end