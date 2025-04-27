function [equalfan, dfan, idealphi, midU_phi] = idealfanangles(fanangles, mid_U, delta_d)
% ideal equal fan and ideal equal radial phi

% size
[Npixel, Nslice] = size(fanangles);
% I know for flying focal, the Npixel is Npixel*Nfocal

% ideal equal fan angle
fanangles_mean = mean(fanangles, 2);
pixelindex = 1:Npixel;
dfan = ((pixelindex-mid_U)*fanangles_mean) / sum((pixelindex-mid_U).^2);
equalfan = dfan.*((1:Npixel)'-mid_U);

% ideal phi (equal radial)
t0 = mod(mid_U, 1);
delta_phi = delta_d*sign(dfan);

t1 = ceil(sin(equalfan(1))/delta_phi + t0);
t2 = floor(sin(equalfan(end))/delta_phi + t0);
idealphi = asin(((t1:t2)'-t0).*delta_phi);
% I know the delta_d = hx_ISO/Nfocal/SID;

% midU of phi
midU_phi = -t1+1+t0;

end