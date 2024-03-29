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
Ny = imgsize(1);
Nx = imgsize(end);
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

% interp target
[msXa, msYa] = meshgrid(Xa, Ya);
Pa = atan2(msYa, msXa)./pi;
Pa = mod(Pa(:).*Ntheta - Za.*flag_rotdirect, Ntheta*2) - Ntheta;

Ra = sqrt(msYa.^2 + msXa.^2);
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
Va = reshape(Va, Ny, Nx, Nimg);
Ra = reshape(Ra, Ny, Nx, Nimg);
Pa = reshape(Pa, Ny, Nx, Nimg);
% Ra = repmat(Ra, 1, 1, Nimg);
% Pa = repmat(Pa, 1, 1, Nimg);

% interp 3D
img = interp3(Vb, Rb, Pb, raw, Va, Ra, Pa);
% done

end
