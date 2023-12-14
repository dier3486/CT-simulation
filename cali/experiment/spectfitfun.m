function [r, resp] = spectfitfun(x, Dp, Dair, spectrange, samplekeV, mu_ref, Dtarg, lambda)

% x = normr(x);
% xt = [0; abs(x(:))]; 
% xt = [0; x(:)]; 
xt = x(:);
Nx = length(xt);
t = linspace(spectrange(1), spectrange(2), Nx);
% cs = spline(t, xt);
% r2 = ppval(cs, C.samplekeV(:));
resp = pchip(t, xt, samplekeV(:));
resp(samplekeV(:)<spectrange(1)) = 0;

% GOS model
% r2 = r2.*C.GOS;

r = (-log(Dp*resp)+log(Dair*resp))./mu_ref - Dtarg;

if nargin > 7
    rr = resp(2:end-1).*2 - resp(1:end-2) - resp(3:end);
    r = [r; sum(rr.^2)./sum(resp.^2).*length(resp).*lambda];
end

end