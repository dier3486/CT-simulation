function [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status)
% recon node, anti ring on image
% [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);

% parameters from pipe
ARprm = prmflow.pipe.(status.nodename);
if isfield(ARprm, 'alpha')
    alpha = ARprm.alpha;
else
    alpha = 1.0;
end
if isfield(ARprm, 'Ntheta')
    Ntheta = ARprm.Ntheta;
else
    Ntheta = 192;
end

% parameters from recon
FOV = prmflow.recon.FOV;
imagesize = prmflow.recon.imagesize;
imagecenter = prmflow.recon.imagecenter;
if isfield(prmflow.recon, 'windowcenter') && isfield(prmflow.recon, 'windowwidth')
    Lb =  prmflow.recon.windowcenter - prmflow.recon.windowwidth/2 + 1000;
    Ub =  prmflow.recon.windowcenter + prmflow.recon.windowwidth/2 + 1000;
else
    Lb = 950;
    Ub = 1150;
end
d_radius = prmflow.recon.delta_d*imagesize/FOV;

% cneter on image space
centerfix = imagecenter(:, 1:2).*(imagesize/FOV);
% anti ring on image space
imgfix = antiringonimage(dataflow.image, centerfix, Lb, Ub, Ntheta, d_radius);
% add to image
dataflow.image = dataflow.image - imgfix.*alpha;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end