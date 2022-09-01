function ImgOut = getBoneImg(ImgIn, BoneCurve)

% % old codes
% minValue = min(BoneCurve(:,1));
% maxValue = max(BoneCurve(:,1));
% ImgIn(ImgIn < minValue) = 0;
% ImgIn(ImgIn > maxValue) = maxValue;
% idx = find(ImgIn > 0);
% ImgIn(idx)=interp1(BoneCurve(:,1), BoneCurve(:,2), ImgIn(idx));
% ImgOut = ImgIn;

minValue = min(BoneCurve(:,1));
ImgOut = interp1(BoneCurve(:,1), BoneCurve(:,2), ImgIn, 'linear', 'extrap');
ImgOut = ImgOut.*(ImgIn>minValue);

end
