function imgfix = antiringonimage(img, center, Lb, Ub)
% anti-ring on image space
% img = antiringonimage(img, center);

% gray cut
img(img<Lb) = Lb;
img(img>Ub) = Ub;

% r-theta
raw = rthetatrans(img, center);
Nth = size(raw, 2);

% rawfix
rawfix = raw(2:end-1, :, :) - raw(1:end-2, :, :)./2 - raw(3:end, :, :)./2;
rawfix = mean(rawfix, 2, 'omitnan');
rawfix = rawfix + flipud(rawfix);

% radius cut
Na = size(img, 1);
Nb = size(rawfix, 1);
Ncut = max(ceil((Nb-Na)/2), 0);
rawfix(1:Ncut, :, :) = 0;
rawfix(end-Ncut+1:end, :, :) = 0;

% rep
rawfix = repmat(rawfix, 1, Nth);

% inv r-theta
imgfix = rthetainv(rawfix, center);

end