function t = axialchordhomeomorph(eta, zeta, Nslice)
% prepare for homeomorphic axial reconstruction of 'chord' beam 
% t = axialchordhomeomorph(eta, zeta, Nslice);

% inputs
dsize = size(eta);
if dsize(end)==1
    dsize = dsize(1:end-1);
end
ZgridA = -Nslice+1:Nslice;
eta = eta(:);
zeta = zeta(:);

% dz and gap
% 0
dz_0 = (1+zeta./sqrt(1-eta.^2));
% pi
dz_pi = (1-zeta./sqrt(1-eta.^2));
% gap
gap = ones(size(eta));

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