function rebin = rebinprepare(detector, focalposition, Nviewprot, isQDO)
% rebin prepare
% recon = rebinprepare(detector, focalposition, Nview, isQDO);
% where the inputs,
%   detector is the struct of detector corr, e.g. prmflow.system.detector,
%   focalposition is the focal position, 
%     e.g. prmflow.system.focalposition(ifocal, :),
%     fly-focal is not supported yet,
%   Nviewprot is the view number per rotation, e.g. prmflow.recon.Nviewprot,
%   isQDO it the bool value of is QDO.
% The returns are,
%   rebin.delta_view,         delta view angle
%   rebin.interalpha_azi,     interp coeeficients for Azi-rebin
%   rebin.vindex1_azi,
%   rebin.vindex2_azi,        index for Azi-rebin
%   rebin.Npixel,             Npixel or Npixel*2 for QDO
%   rebin.Nviewprot,          Nviewprot or Nviewprot/2 for QDO
%   rebin.QDOorder,           the alternating order of the pixels for QDO
%   rebin.interalpha_rad,     interp coeeficients radial-rebin
%   rebin.radialindex,        index for radial-rebin
%   rebin.Nreb,               pixel number after rebin
%   rebin.delta_d,            pixle size after rebin
%   rebin.midchannel,         midchannel after rebin

% default, not QDO
if nargin<4
    isQDO = false;
end

% parameters to use
Npixel = double(detector.Npixel);
Nslice = double(detector.Nmergedslice);
mid_U = single(detector.mid_U);
Nps = Npixel*Nslice;
hx_ISO = detector.hx_ISO;

% fan angles
y = detector.position(1:Npixel, 2) - focalposition(2);
x = detector.position(1:Npixel, 1) - focalposition(1);
fanangles = atan2(y, x);
% I know the fanangles of each slice are equal

% perpare for Azi rebin
delta_view = pi*2/Nviewprot;
f = fanangles./delta_view;
viewindex = double(floor(f));
rebin.delta_view = delta_view;
rebin.interalpha_azi = repmat(f-viewindex, Nslice, 1);
viewindex = viewindex + 1;  % start from 0
rebin.startvindex = mod(max(viewindex), Nviewprot)+1;
viewindex = repmat(viewindex, Nslice, Nviewprot) + repmat(0:Nviewprot-1, Nps, 1);
rebin.vindex1_azi = mod(viewindex-1, Nviewprot).*Nps + repmat((1:Nps)', 1, Nviewprot);
rebin.vindex2_azi = mod(viewindex, Nviewprot).*Nps + repmat((1:Nps)', 1, Nviewprot);

% prepar for radial rebin
if isQDO
    % QDO order
    [a1, a2] = QDOorder(Npixel, mid_U);
    s1 = ~isnan(a1);
    s2 = ~isnan(a2);
    rebin.Npixel = max([a1, a2]);
    rebin.QDOorder = [a1(:), a2(:)];
    % d0 is the distance from ray to ISO
    d0 = -detector.SID.*cos(fanangles);
    % QDO d
    d = nan(rebin.Npixel, 1);
    d(a1(s1)) = d0(s1);
    d(a2(s2)) = -d0(s2);
    % delta_t and mid_t
    delta_d = hx_ISO/2.0;
    mid_t = 0.5;
    % Nview
    rebin.Nviewprot = Nviewprot/2;
else
    rebin.Npixel = Npixel;
    d = -detector.SID.*cos(fanangles);
    delta_d = detector.hx_ISO;
    mid_t = mod(detector.mid_U, 1);
    % Nview
    rebin.Nviewprot = Nviewprot;
end
% radial interp grid
t1 = ceil(min(d)/delta_d + mid_t);
t2 = floor(max(d)/delta_d + mid_t);
Nreb = t2-t1+1;
tt = ((t1:t2)-mid_t)'.*delta_d;
% interp index
fd = d./delta_d + mid_t;
dindex = floor(fd) - t1 + 2;
dindex(dindex<=0) = 1;
dindex(dindex>Nreb) = Nreb+1;
tindex = nan(Nreb+1, 1);
tindex(dindex) = 1:rebin.Npixel;
tindex = fillmissing(tindex(1:end-1), 'previous');
% got it
rebin.Nreb = Nreb;
rebin.delta_d = delta_d;
rebin.radialindex = tindex;
rebin.interalpha_rad = (tt - d(tindex))./(d(tindex+1)-d(tindex));
% midchannel
rebin.midchannel = -t1+1+mid_t;

end