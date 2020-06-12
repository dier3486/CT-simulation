function p = pinprojectfun(viewangle, Nslice, hz, varargin)
% pin projection model function
% p = pinprojection(viewangle, Nslice, phi, r, phi, zeta_x, zeta_y, hz, rotscale, midc, dslope, dscale, isooff, isophase, zshift);
% r is r/SID, zeta is zeta/SID.

Nvin = 16;
if nargin<Nvin
    varargin{Nvin - 3} = 0;
    varargin(max(nargin - 2, 1) :end) = {0};
end

[r, phi, zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, isooff2, isophase2, zshift] = varargin{:};

% skip zshift, fix in 0
zshift = 0;

% test to fix something
% dslope = 0;
% zeta_x = 0; zeta_y = 0;

zeta_x = zeta_x/1000;
zeta_y = zeta_y/1000;

% fix view angle
viewangle = viewangle+[0 cumsum(diff(viewangle)<0)].*(pi*2);
viewangle = (viewangle - viewangle(1)).*(rotscale+1) + viewangle(1);

% view woble
% viewangle = viewangle + sin((viewangle+viewphase)*2).*viewoff;

% fix mid
S = (1:Nslice)' - (Nslice+1)/2 + zshift;
midfix = midc + S.*dslope + sin(viewangle+isophase).*isooff + sin((viewangle+isophase2).*2).*isooff2;
% midfix = midc + S.*dslope + sin(viewangle+isophase).*isooff + sin((viewangle+isophase + pi/8).*2).*isooff2;
% midfix = midc + S.*dslope;

p = projectonslices(viewangle, r, phi, Nslice, zeta_x, zeta_y, zshift, hz).*(dscale+1) + midfix;



end

function p = projectonslices(viewangle, r, phi, Nslice, zeta_x, zeta_y, zshift, hz)
% pin projection model function
% p = pinprojection(viewangle, r, phi, Nslice, zeta_x, zeta_y, hz);
% r is r/SID, zeta is zeta/SID.

% S should be [... -3/2 -1/2 1/2 3/2 ...] to employ the z-steps of the slices
S = (1:Nslice)' - (Nslice+1)/2 + zshift;
% L is S * D where D is the length from focal spot to pin center of each view
L = S * sqrt(1 + r^2 + sin(phi-viewangle).*(r*2));

% (phi, r) is the polar coorinates of the pin position (x, y)
% but for each slice the pin position (x, y) moved to (x+delta_x, y+delta_y) due to the pin is slope 
% they read,
phi_s = atan2(r*sin(phi) + L.*zeta_y.*hz, r*cos(phi) + L.*zeta_x.*hz);
r_s = sqrt((r*cos(phi) + L.*zeta_x.*hz).^2 + (r*sin(phi) + L.*zeta_y.*hz).^2);

% pin projection angle
p = atan2(-cos(phi_s-viewangle), (1./r_s+sin(phi_s-viewangle)));

end