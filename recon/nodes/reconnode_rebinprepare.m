function [prmflow, status] = reconnode_rebinprepare(prmflow, status)
% recon node, rebin prepare
% [prmflow, status] = reconnode_rebinprepare(prmflow, status);

% parameters to use in prmflow
Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.recon.Nfocal;
rebinpipe = prmflow.pipe.(status.nodename);

% isQDO
if isfield(rebinpipe, 'QDO')
    isQDO = rebinpipe.QDO;
else
    isQDO = false;
end
% debug
% isQDO = true;

% detector
detector = prmflow.system.detector;

% fan angles & focal angle(s)
if isfield(prmflow.recon, 'fanangles')
    fanangles = prmflow.recon.fanangles;
    focalangle = prmflow.recon.focalangle;
else
    [fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
    prmflow.recon.fanangles = fanangles;
    prmflow.recon.focalangle = focalangle;
end

% rebin prepare
prmflow.rebin = rebinprepare(detector, fanangles, focalangle, Nviewprot, isQDO);
prmflow.rebin.isQDO = isQDO;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end