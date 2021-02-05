function t = axialconehomeomorph(eta, zeta, Nslice, theta_tilt)
% prepare for homeomorphic axial reconstruction of conebeam 
% t = axialconehomeomorph(eta, zeta, Nslice, theta_tilt);

% inputs
if nargin<4
    theta_tilt = 0;
end
dsize = size(eta);
if dsize(end)==1
    dsize = dsize(1:end-1);
end
ZgridA = -Nslice+1:Nslice;
eta = eta(:);
zeta = zeta(:);

% dz and gap
s_tilt = abs(eta)<abs(sin(theta_tilt));
% 0
dz_0 = zeros(size(eta));
dz_0(s_tilt) = (1+zeta(s_tilt)./sqrt(1-eta(s_tilt).^2));
dz_0(~s_tilt) = (sqrt(1-eta(~s_tilt).^2)+zeta(~s_tilt))./cos(theta_tilt);
% pi
dz_pi = zeros(size(eta));
dz_pi(s_tilt) = (1-zeta(s_tilt)./sqrt(1-eta(s_tilt).^2));
dz_pi(~s_tilt) = (sqrt(1-eta(~s_tilt).^2)-zeta(~s_tilt))./cos(theta_tilt);
% gap
gap = ones(size(eta));
gap(~s_tilt) = Nslice - sqrt(1-eta(~s_tilt).^2).*((Nslice-1)/cos(theta_tilt));

% t
t = nan([dsize Nslice*2], class(eta));
% 0
s_0 = dz_0 >= abs(ZgridA.*2-1)./(Nslice-1);
t_tmp = (ZgridA-1/2)./dz_0+1/2;
t(s_0) = t_tmp(s_0);
% pi
s_pi = dz_pi >= (Nslice*2-abs(ZgridA.*2-1))./(Nslice-1);
t_tmp = Nslice.*sign(ZgridA-1/2)+1/2 + (-Nslice.*sign(ZgridA-1/2)+ZgridA-1/2)./dz_pi;
t(s_pi) = t_tmp(s_pi);
% gap
s_gap = ~s_0 & ~s_pi;
t_tmp = ((Nslice-1)/2).*sign(ZgridA-1/2)+1/2 + (ZgridA-1/2 - dz_0.*(((Nslice-1)/2).*sign(ZgridA-1/2)))./gap;
t(s_gap) = t_tmp(s_gap);

t = reshape(t, [dsize Nslice*2]);
end