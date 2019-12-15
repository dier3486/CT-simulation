function [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status)
% recon node, housefield correction
% [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status);

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.protocol.viewperrot;
focalspot = prmflow.system.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.system.Nfocal;
% fly-focal is not supported yet

% only Axial, one shot

% detector
detector = prmflow.system.detector;
mid_U = detector_corr.mid_U;
hx_ISO = detector_corr.hx_ISO;
Nps = Npixel*Nslice;

% fan angles
y = detector.position(1:Npixel, 2) - focalposition(2);
x = detector.position(1:Npixel, 1) - focalposition(1);
fanangles = atan2(y, x);
% I know the fanangles of each slice are equal

% rebin 1
delta_view = pi*2/Nviewprot;
f = fanangles./delta_view;
viewindex = double(floor(f));
interalpha = repmat(f-viewindex, Nslice, 1);
viewindex = viewindex + 1;  % start from 0
startvindex = mod(max(viewindex), Nview)+1;
viewindex = repmat(viewindex, Nslice, Nview) + repmat(0:Nview-1, Nps, 1);
vindex1 = mod(viewindex-1, Nview).*Nps + repmat((1:Nps)', 1, Nview);
vindex2 = mod(viewindex, Nview).*Nps + repmat((1:Nps)', 1, Nview);

A = zeros(Nps, Nview*Nshot, 'single');
for ishot = 1:Nshot
    start_ishot = (ishot-1)*Nps*Nview;
    viewindex = (1:Nview) + (ishot-1)*Nview;
    A(start_ishot+vindex1) = dataflow.rawdata(:, viewindex).*repmat(1-interalpha, 1, Nview);
    A(start_ishot+vindex2) = A(start_ishot+vindex2) + dataflow.rawdata(:, viewindex).*repmat(interalpha, 1, Nview);
    % start angle for first rebin view
    A(:, viewindex) = [A(:, (startvindex:Nview)+(ishot-1)*Nview) A(:, (1:startvindex-1)+(ishot-1)*Nview)];
end
% start angle
viewangle = dataflow.rawhead.viewangle;
viewangle = [viewangle(startvindex:end) viewangle(1:startvindex-1)];
startviewangle = viewangle(1);

% to return
dataflow.rawdata = A;
dataflow.rawhead.viewangle = viewangle;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end