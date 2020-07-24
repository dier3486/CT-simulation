function img = helicalrthetainv(raw, Nrow, varargin)
% interse helical r-theta transform of the raw
% img = helicalrthetainv(raw, Nrow, center, imgsize, Ntheta, d, flag_rotdirect);
% or, raw = rthetatrans(img, Nrow);

% raw size
[Nb, Ntheta_all, Nslice] = size(raw);

% default inputs
%               center        imgsize   Ntheta   d        flag_rotdirect
defaultinput = {zeros(1, 2),  512,      192,     1.0,     1              };
% input coeffients
[center, imgsize, Ntheta, d, flag_rotdirect] = cleaninputarg(defaultinput, varargin{:});

if floor(Ntheta/Nslice) ~= Ntheta/Nslice
    error('The Ntheta must by integer times to size(raw,3)!');
end

% Nimg
Nimg = round(Ntheta_all/Ntheta/2*Nrow);

% fill up raw
raw = cat(3, raw, flipud(raw(:,:,1)));
if flag_rotdirect>=0
    raw = cat(2, raw, zeros(Nb, Ntheta, Nslice+1));
else
    raw = cat(2, zeros(Nb, Ntheta, Nslice+1), raw);
end
for ii = 1:Nslice
    raw(:, :, ii+1) = circshift(raw(:, :, ii+1), ii*Ntheta/Nslice*flag_rotdirect, 2);
end

% XY grid of image(s)
Nx = imgsize(1);
Ny = imgsize(end);
Xa = (-(Nx-1)/2 : (Nx-1)/2) - center(1, 1);
Ya = (-(Ny-1)/2 : (Ny-1)/2) - center(1, 2);
% NOTE: no movable image center (tilt) in helical
Za = ((1:Nimg) - 1/2).*(Ntheta*2/Nrow);

% theta grid of raw
Vb = 0 : (Ntheta_all -1 + Ntheta);
% radius grid of raw
Rb = (-(Nb-1)/2 : (Nb-1)/2).*d;
% slice grid of raw
Pb = (0:Nslice).*(1/Nslice);

% interp target (ndgrid)
[ndXa, ndYa] = ndgrid(Xa, Ya);
Pa = atan2(ndYa, ndXa)./pi;
Pa = mod(Pa(:).*Ntheta - Za.*flag_rotdirect, Ntheta*2) - Ntheta;

Ra = sqrt(ndYa.^2 + ndXa.^2);
Ra = Ra(:).*sign(-Pa);

Pa(Pa<0) = Pa(Pa<0) + Ntheta;
% Va
if flag_rotdirect>0
    Va = Za + Pa;
elseif flag_rotdirect<0
    Va = Za + (Ntheta-Pa);
end
% norm Pa to [0 1]
Pa = Pa./Ntheta;

% repmat and reshape
Va = reshape(Va, Nx, Ny, Nimg);
Ra = reshape(Ra, Nx, Ny, Nimg);
Pa = reshape(Pa, Nx, Ny, Nimg);
% Ra = repmat(Ra, 1, 1, Nimg);
% Pa = repmat(Pa, 1, 1, Nimg);

% interp 3D
img = interp3(Vb, Rb, Pb, raw, Va, Ra, Pa);
% done

end
