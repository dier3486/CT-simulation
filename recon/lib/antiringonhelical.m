function imgfix = antiringonhelical(img, center, Nrow, rotdirect, Lb, Ub, varargin)
% anti-ring on image space
% imgfix = antiringonhelical(img, center, Nrow, rotdirect, Lb, Ub, Ntheta, d, Zsample, flag_even, restcut, ringfilter, Nsect, sectmethod);
% or imgfix = antiringonhelical(img, center, Nrow, rotdirect, Lb, Ub);
% then the fixed img = img-imgfix;
% INPUT:
%   img             the images to fix, e.g. a matrix in size 512x512xNimg (x,y,z)
%   center          the XY position of ISO center on images space
%   Nrow            the image number in one pitch
%   rotdirect       the rotation direction, 1 or -1,  1 means x->y, -1 means y->x
%   Lb, Ub          the window [Lb Ub] to remove rings
%   ------ 
%   Ntheta          the 'view' nubmer per half roation in r-theta ransform, default value is 192
%   d               the sample step length in radius, default value is 1
%   Zsample         the sample number in Z direction, Nrow*Zsample is the 'slice' number, usually Zsample>=1/pitch.
%   flag_even       the flag of the symmetry of radius samples, true is even-symmetry, false is odd-symmetry, default value true
%   restcut         coeffient to control the restaint of ring fix, in 0 to 1, 0 is all pass, 1 is all block, default value is 0.1
%   ringfilter      the filter to pick up rings, default is [-1/2 1 -1/2]
%   Nsect           the section number of anti-ring per roation
%   sectmethod      the method to link the rings between the sections, 'linear' or 'spline', default value is 'spline'.

% default inputs
%               Ntheta, d,      Zsample,    flag_even,  restcut,    ringfilter,     Nsect,  sectmethod
defaultinput = {192,    1.0,    2,          true,       0.1,        [-1/2 1 -1/2],  4,      'spline'};
% input coeffients
[Ntheta, d, Zsample, flag_even, restcut, ringfilter, Nsect, sectmethod] = cleaninputarg(defaultinput, varargin{:});

% gray cut
img(img<Lb) = Lb;
img(img>Ub) = Ub;

% r-theta
raw = helicalrthetatrans(img, Nrow, center, Ntheta, d, Zsample, flag_even, rotdirect);
[Nb, Nv, Nslice] = size(raw);

% radius cut
imagesize = size(img, [1 2]);
Na = max(imagesize);
% Na = size(img, 1);
Ncut = max(ceil( (Nb - Na/d)/2 ), 0);

% ring filter
raw = reshape(raw, Nb, []);
rawring = conv2(raw, ringfilter(:), 'same');

% restaint
fixrest = [zeros(1, size(raw, 2)); diff(raw)];
fixrest = abs(fixrest) + flipud(abs(fixrest));
fixrest = 1 - (fixrest./(Ub-Lb).*restcut).^2;
fixrest(fixrest<0) = 0;
% apply
rawring = rawring.*fixrest;

% full
rawring = reshape(rawring, Nb, Nv, Nslice);

% loop sections to 
Nvpersect = floor(Ntheta*2/Nsect);
Nsect_all = ceil(Nv/Nvpersect);
ringact = zeros(Nb, Nsect_all+1, Nslice);
for isect = 1:Nsect_all
    index = (-Nvpersect : Nvpersect) + Nvpersect*(isect-1) + 1;
    index = index((index>0) & (index<=Nv));
    ringact(:, isect, :) = median(rawring(:, index, :), 2, 'omitnan');
end
ringact(:, Nsect_all+1, :) = ringact(:, Nsect_all, :);

% interpolation
view = 0 : Nv - 1;
viewact = (0:Nsect_all).*Nvpersect;
rawfix = zeros(Nb, Nv, Nslice);
active_index = Ncut+1 : Nb-Ncut;
switch sectmethod
    case 'linear'
    % linear interp
        [intp_index, intp_alpha] = interpprepare(viewact, view, 'extrap');
        rawfix(active_index, :, :) = ringact(active_index, intp_index(:, 1), :).*intp_alpha(:, 1)' + ...
            ringact(active_index, intp_index(:, 2), :).*intp_alpha(:, 2)';
    case 'spline'
        % spline interp
        for irow = 1 : Nslice
            for ii = active_index
                ringact_ii = [0 ringact(ii, :, irow) 0];
                rawfix(ii, :, irow) = spline(viewact, ringact_ii, view);
            end  
        end
    otherwise
        error('Unknown interp method %s!', sectmethod);
end
rawfix = rawfix.*2;

% inv r-theta
imgfix = helicalrthetainv(rawfix, Nrow, center, imagesize, Ntheta, d, rotdirect);

end