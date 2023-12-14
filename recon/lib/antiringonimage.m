function imgfix = antiringonimage(img, center, Lb, Ub, varargin)
% anti-ring on image space
% img = antiringonimage(img, center, Lb, Ub, Ntheta, d, flag_even, restcut, ringfilter, Nsect, sectmethod);
% or, img = antiringonimage(img, center, Lb, Ub);

% gray cut
img(img<Lb) = Lb;
img(img>Ub) = Ub;

% default inputs
%               Ntheta, d,      flag_even,  restcut,    ringfilter,     Nsect,  sectmethod
defaultinput = {192,    1.0,    true,       0.1,        [-1/2 1 -1/2],  4,      'spline'};
% input coeffients
[Ntheta, d, flag_even, restcut, ringfilter, Nsect, sectmethod] = cleaninputarg(defaultinput, varargin{:});

% r-theta
raw = rthetatrans(img, center, Ntheta, d, flag_even);
Nb = size(raw, 1);
% Nth = size(raw, 2);
raw = reshape(raw, Nb, []);

% radius cut
% imagesize = size(img, [1 2]); % matlab 2018+
imagesize = [size(img, 1) size(img, 2)];
Na = max(imagesize);
% Nb = size(rawfix, 1);
Ncut = max(ceil((Nb-Na/d)/2), 0);

% ring filter
rawring = conv2(raw, ringfilter(:), 'same');

% restaint
fixrest = [zeros(1, size(raw, 2), 'like', img); diff(raw)];
fixrest = abs(fixrest) + flipud(abs(fixrest));
fixrest = 1 - (fixrest./(Ub-Lb).*restcut).^2;
fixrest(fixrest<0) = 0;
% apply
rawring = rawring.*fixrest;

% full
rawring = reshape(rawring, Nb, Ntheta, []);
rawring = [rawring flipud(rawring)];

% loop sections to 
Nimg = size(img, 3);
ringact = zeros(Nb, Nsect+1, Nimg, 'like', img);
for isect = 1:Nsect
    Nv = floor(Ntheta*2/Nsect);
    index = (1 : Nv*2+1) + Nv*(isect-1);
    index = mod(index-1, Ntheta*2)+1;
    ringact(:, isect+1, :) = median(rawring(:, index, :), 2, 'omitnan');
end
ringact(:, 1, :) = ringact(:, Nsect+1, :);

% interpolation
theta = cast((0:Ntheta-1).*(pi/Ntheta), 'like', img);
thetaact = cast(linspace(0, pi*2, Nsect+1), 'like', img);
rawfix = zeros(Nb, Ntheta, Nimg, 'like', img);
active_index = Ncut+1:Nb-Ncut;
switch sectmethod
    case 'linear'
    % linear interp
        [intp_index, intp_alpha] = interpprepare(thetaact, theta, 'extrap');
        rawfix(active_index, :, :) = ringact(active_index, intp_index(:, 1), :).*intp_alpha(:, 1)' + ...
            ringact(active_index, intp_index(:, 2), :).*intp_alpha(:, 2)';
    case 'spline'
        % spline interp
        for irow = 1:Nimg
            for ii = active_index
                ringact_ii = [0 ringact(ii, :, irow) 0];
                rawfix(ii, :, irow) = spline(thetaact, ringact_ii, theta);
            end  
        end
    otherwise
        error('Unknown interp method %s!', sectmethod);
end
rawfix = rawfix.*2;

% inv r-theta
imgfix = rthetainv(rawfix, center, imagesize, d);

end