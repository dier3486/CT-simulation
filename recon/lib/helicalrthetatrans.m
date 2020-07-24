function raw = helicalrthetatrans(img, Nrow, varargin)
% helical r-theta transform of the img
% raw = helicalrthetatrans(img, Nrow, center, Ntheta, d, Zsample, flag_even, flag_rotdirect);
% or, raw = helicalrthetatrans(img, Nrow);

% image size
[Nx, Ny, Nimg] = size(img);

% default inputs
%               center        Ntheta   d        Zsample   flag_even   flag_rotdirect
defaultinput = {zeros(1, 2),  192,     1.0,     2,        true,       1              };
% input coeffients
[center, Ntheta, d, Zsample, flag_even, flag_rotdirect] = cleaninputarg(defaultinput, varargin{:});
% only integer Zample
Zsample = floor(Zsample);
% slice number
Nslice = Nrow*Zsample;
if floor(Ntheta/Nslice) ~= Ntheta/Nslice
    % not common?
    Ntheta = ceil(Ntheta/Nslice)*Nslice;
    % come on common nmnmnmm
    warning('The Ntheta is changed to %d to common with the Nrow*Zsample!', Ntheta);
end

% XY grid of image
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(1, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(1, 2);
% NOTE: no movable image center (tilt) in helical
% Z grid of image (add boundry)
Za = 0:Nimg+1;

% theta grid of raw
Ntheta_all = floor(2*Ntheta*Nimg/Nrow);
% I know the Ntheta the number of 'views' per half rotation
Vb = (0 : Ntheta_all-1).*(pi/Ntheta).*flag_rotdirect;
% radius grid of raw
Nb = ceil(max(Nx, Ny)/d*sqrt(2)/2)*2;
if ~flag_even
	Nb = Nb+1;
end
Rb = (-(Nb-1)/2 : (Nb-1)/2).*d;
% Z grid to shift Vb

% Zgrid = 0 : 1/Zsample : Nrow-1/Zsample;
% Vbshift = Zgrid.*(2*pi/Nrow).*flag_rotdirect;
Vbshift = (0 : Nslice-1).*(pi/Nslice);
Vb = Vb(:) + Vbshift;

% interp target (ndgrid)
Xb = Rb(:) * cos(Vb(:)');
Yb = Rb(:) * sin(Vb(:)');
Zb = repmat((0 : Ntheta_all-1).*(Nrow/Ntheta/2) + 1/2, Nb, Nslice);
% reshape
Xb = reshape(Xb, Nb, Ntheta_all, Nslice);
Yb = reshape(Yb, Nb, Ntheta_all, Nslice);
Zb = reshape(Zb, Nb, Ntheta_all, Nslice);

% img boundry
img = cat(3, img(:,:,1).*2-img(:,:,2), img, img(:,:,end).*2-img(:,:,end-1));
% interp 3D
raw = interp3(Xa, Ya, Za, permute(img, [2 1 3]), Xb, Yb, Zb);
% done

end
