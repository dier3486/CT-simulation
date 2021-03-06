function [prmflow, status] = reconnode_rebinprepare(prmflow, status)
% recon node, rebin prepare
% [prmflow, status] = reconnode_rebinprepare(prmflow, status);

% parameters to use in prmflow
Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.recon.Nfocal;
% fly-focal is not supported yet
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

% rebin prepare
prmflow.rebin = rebinprepare_bk0604(detector, focalposition, Nviewprot, isQDO);
prmflow.rebin.isQDO = isQDO;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end