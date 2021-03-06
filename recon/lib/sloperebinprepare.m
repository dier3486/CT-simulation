function rebin = sloperebinprepare(detector, fanangles, focalangle, Nviewprot, gantrytilt)
% 'slope' rebin prepare for Axial
% support QDO, XDFS and QDO+XDFS, but not support Z-DFS.
% recon = rebinprepare(detector, focalposition, Nview, isQDO);
% where the inputs,
%   detector,                   the struct of detector corr, e.g. prmflow.system.detector,
%   fanangles,                  they are,
%   focalangle,                 [fanangles, focalangle] = detpos2fanangles(detposition, focalposition);                              
%   Nviewprot,                  the view number per rotation, e.g. prmflow.recon.Nviewprot,
%   gantrytilt,                 gantry tilt
% The returns are,
%   rebin.delta_view,           delta view angle
%   rebin.faninterpkern,        interp coeeficients of equal fan-angles for slope fan-Radial rebin
%   rebin.dfan,                 equal fan size
%   rebin.idealphi,             equal radial fan-ganles (ideal phi)
%   rebin.Npixel,               Npixel
%   rebin.Nviewprot,            Nviewprot after rebin
%   rebin.Nreb,                 pixel number after rebin
%   rebin.delta_d,              pixle size after rebin
%   rebin.midchannel,           midchannel after rebin
%   rebin.midU,                 midchannel before slope fan-Radial rebin
%   rebin.Yshift,               equivalent Y shift in tilt
%   rebin.Zgrid,                interp coeeficients of equivalent tilt for slope fan-Radial rebin
%   rebin.DFSviewinterp,        DFS

% prepare subfun (Inputs: detector, fanangles, focalangle, Nviewprot, gantrytilt)
[Npixel, Nslice, ~] = size(fanangles);
Nps = Npixel*Nslice;
midU = detector.mid_U;
hx_ISO = detector.hx_ISO;
SID = detector.SID;
Nfocal = length(focalangle);
Npixelmf = Npixel*Nfocal;
% fanangles-focalangle
% fanangles = reshape(reshape(fanangles, Nps, []) - focalangle, Npixel, Nslice, Nfocal);
fanangles = reshape(fanangles, Nps, []) - focalangle;
fanangles = reshape(fanangles', Npixelmf, Nslice);

% mid_U
if Nfocal == 2
    midU = DFSmidchannel(midU(1), abs(focalangle(1)-pi/2) > abs(focalangle(2)-pi/2));
else
    midU = midU(1);
end
% delta d
delta_d = hx_ISO/Nfocal;

% ideal equal fan angle and ideal equal radial phi
[equalfan, dfan, idealphi, midU_phi] = idealfanangles(fanangles, midU, delta_d/SID);
Nreb = length(idealphi);

faninterpkern = zeros(Npixelmf, Nslice, 'single');
for islice = 1:Nslice
    faninterpkern(:, islice) = interp1(fanangles(:,islice), 1:Npixelmf, equalfan, 'linear', 'extrap');
end

% Yshift by tilt
Zdet = mean(reshape(detector.position(:,3), Npixel, Nslice));
Yshift = -Zdet.*tan(gantrytilt)./detector.SDD;

% Zgrid
Zsec = cos(gantrytilt)./cos(idealphi);
Zsec(Zsec>1) = 1;
Zgrid = Zsec*(-(Nslice-1)/2 : (Nslice-1)/2) + (Nslice+1)/2;

% DFS
if Nfocal == 2
    dview = pi*2/Nviewprot;
    DFSviewinterp = -((focalangle-pi/2)./dview + [-1/2 1/2])./2;
else
    DFSviewinterp = [];
end

% prepare for fan-Radial
rebin.Nviewprot = Nviewprot/Nfocal;
rebin.Nreb = Nreb;
rebin.delta_d = delta_d;
rebin.midchannel = midU;
rebin.midU_phi = midU_phi;
rebin.faninterpkern = faninterpkern;
rebin.dfan = dfan;
rebin.idealphi = idealphi;
rebin.Yshift = Yshift;
rebin.Zgrid = Zgrid;
rebin.DFSviewinterp = DFSviewinterp;

% other prms from detector
rebin.SID = detector.SID;
rebin.delta_z = detector.hz_ISO;


end


