function y = slidecurve(ta, tb, sigma)
% a slide curve
% ta = x-a; tb = b-x;
% 0 free

t = fillmissing((ta-tb)./(ta+tb), 'constant', 0);
y = (ta-tb)./2 + (ta+tb)./2.*sign(t).*sigma.*log((1+exp(1/sigma).*exp(-abs(t)./sigma))./(exp(1/sigma)+exp(-abs(t)./sigma)));

end
