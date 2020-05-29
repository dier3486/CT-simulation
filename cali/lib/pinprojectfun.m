function p = pinprojectfun(viewangle, r, phi, Nslice, zeta_x, zeta_y, hz)
% pin projection model function
% p = pinprojection(viewangle, r, phi, Nslice, zeta_x, zeta_y, hz);
% r is r/SID, zeta is zeta/SID.

if nargin<4
    Nslice = 1;
end

if nargin<5
    zeta_x = 0;
    zeta_y = 0;
    hz = 1.0;
end

% test
% SID = 570;
% hz = 0.55;
% 
% r0 = 200;
% phi0 = pi/4 + 0.1;
% 
% r1 = r0/SID;
% theta = linspace(0, pi*2, 200);
% 
% pt1 = atan2(-cos(phi0-theta), (1/r1+sin(phi0-theta)));
% 
% zeta = randn(1, 2)./SID;
% Nslice = 16;
% S = (1:Nslice)' - (Nslice+1)/2;
% 
% L1 = sqrt(1 + r1^2 + sin(phi0-theta).*(r1*2));
% phi = atan2(r1*sin(phi0) + (S*L1).*zeta(2).*hz, r1*cos(phi0) + (S*L1).*zeta(1).*hz);
% r = sqrt((r1*cos(phi0) + (S*L1).*zeta(1).*hz).^2 + (r1*sin(phi0) + (S*L1).*zeta(2).*hz).^2);
% 
% pt2 = atan2(-cos(phi-theta), (1./r+sin(phi-theta)));

S = (1:Nslice)' - (Nslice+1)/2;
L = S * sqrt(1 + r^2 + sin(phi-viewangle).*(r*2));
phi_s = atan2(r*sin(phi) + L.*zeta_y.*hz, r*cos(phi) + L.*zeta_x.*hz);
r_s = sqrt((r*cos(phi) + L.*zeta_x.*hz).^2 + (r*sin(phi) + L.*zeta_y.*hz).^2);
p = atan2(-cos(phi_s-viewangle), (1./r_s+sin(phi_s-viewangle)));