function r = spectfitfun(x, Dp, Dair, spectrange, samplekeV, mu_ref, Dtarg)


% x = normr(x);
% xt = [0; abs(x(:))]; 
xt = [0; x(:)]; 
Nx = length(xt);
t = linspace(spectrange(1), spectrange(2), Nx);
% cs = spline(t, xt);
% r2 = ppval(cs, C.samplekeV(:));
r2 = pchip(t, xt, samplekeV(:));
r2(samplekeV(:)<spectrange(1)) = 0;

% GOS model
% r2 = r2.*C.GOS;

r = (-log(Dp*r2)+log(Dair*r2))./mu_ref - Dtarg;
