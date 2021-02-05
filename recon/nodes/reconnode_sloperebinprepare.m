function [prmflow, status] = reconnode_sloperebinprepare(prmflow, status)
% recon node, rebin prepare
% [prmflow, status] = reconnode_sloperebinprepare(prmflow, status);

% parameters to use in prmflow
Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
gantrytilt = prmflow.recon.gantrytilt;
% gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.recon.Nfocal;
rebinpipe = prmflow.pipe.(status.nodename);

if isfield(rebinpipe, 'viewblock') && ~isempty(rebinpipe.viewblock)
    viewblock = rebinpipe.viewblock;
else
    viewblock = 32;
end

% detector
detector = prmflow.system.detector;

% fan angles & focal angle(s)
if isfield(prmflow.recon, 'fanangles')
    fanangles = prmflow.recon.fanangles;
    focalangle = prmflow.recon.focalangle;
else
    [fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
%     prmflow.recon.fanangles = reshape(fanangles, prmflow.recon.Npixel, prmflow.recon.Nslice);
    fanangles = reshape(fanangles, prmflow.recon.Npixel, prmflow.recon.Nslice, []);
    prmflow.recon.fanangles = fanangles;
    prmflow.recon.focalangle = focalangle;
end

% rebin prepare
rebin = sloperebinprepare(detector, fanangles, focalangle, Nviewprot, gantrytilt);
rebin.viewblock = viewblock;

% rebin to prmflow
prmflow.rebin = rebin;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end