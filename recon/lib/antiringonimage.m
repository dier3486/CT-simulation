function ringarc = antiringonimage(img, varargin)
% sparse the ring artifacts of the images, key sub-function of image space anti-ring node.
% ringarc = antiringonimage(img, center, Ntheta, d, Rfov, flag_even, ringcut, ringmover, Nsect);
% or, ringarc = antiringonimage(img);

% default inputs (don't do this in deployments)
%               center,     Crange,     Ntheta, d,      Rfov    flag_even,  ringcut,    ringmover,     Nsect
defaultinput = {[0 0],      [-inf inf], 192,    1.0,    inf,    true,       20.0,       [5 5 5],       1    };
% input coeffients
[center, Crange, Ntheta, d, Rfov, flag_even, ringcut, ringmover, Nsect] = cleaninputarg(defaultinput, varargin{:});

% imagesize
imagesize = size(img, [2 1]);

% gray cut
img(img<Crange(1)) = Crange(1);
img(img>Crange(2)) = Crange(2);

% r-theta
raw = rthetatrans(img, center, Ntheta, d, Rfov, flag_even);
Nb = size(raw, 1);
% I know the Nb <= ceil(Rfov/d)*2 + 1

% movmean on angle direction
m_movm = ringmover(1);
raw = cat(2, flipud(raw(:,end-m_movm+1:end, :)), raw, flipud(raw(:,1:m_movm, :)));
% the raw is in size of (Nb, Ntheta+m_movm*2, Nimage)
raw = movmean(raw, m_movm*2+1, 2);
raw = raw(:, m_movm+1:end-m_movm, :);

% movmedian on radius direction
m_md = ringmover(2);
raw = raw - movmedian(raw, m_md, 1);
% do not use the medfilt1

% cut
raw(raw>ringcut) = ringcut;
raw(raw<-ringcut) = -ringcut;

% full the rawring
Nv = floor(Ntheta/Nsect);
Nv2 = Nv*2 + 1;
Nb_half = ceil(Nb/2);
rawring = cat(2, flipud(raw(end-Nb_half+1:end, end-Nv+1:end, :)), raw(1:Nb_half, :, :), ...
                 flipud(raw(end-Nb_half+1:end, :, :)), raw(1:Nb_half, 1:Nv, :));
% the rawring is in size of (Nb_half, Ntheta*2+Nv*2, Nimage)

% movmean by section
m_acp = ringmover(3);
ringacp = sign(rawring);
rawring = movmean(rawring, Nv2, 2);
ringacp = ringacp.*sign(rawring);
ringacp = ringacp.*(ringacp>0);
ringacp = movmean(ringacp, m_acp*2+1, 2);
rawring = rawring.*ringacp;

% full the rawring (overwriting to raw)
if flag_even
    % even
    raw = cat(1, rawring(:, Nv+1: Nv+Ntheta, :), flipud(rawring(:, Nv+Ntheta+1: Nv+Ntheta*2, :)));
else
    % odd
    raw = cat(1, rawring(:, Nv+1: Nv+Ntheta, :), flipud(rawring(1:end-1, Nv+Ntheta+1: Nv+Ntheta*2, :)));
end
% the rawfix is in size of 

% blunt (hard code)
raw = raw + raw(:, [2:end end], :).*0.5 + raw(:, [1 1:end-1], :).*0.5;

% inv r-theta
ringarc = rthetainv(raw, center, imagesize, d);

end