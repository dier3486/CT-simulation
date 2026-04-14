function cimage = protocol2poster(cimage, protocol)
% the image parameters for postprocesser.
%   prmflow.image = protocol2poster(prmflow.image, prmflow.protocol);
% Note: mostly the prmflow.image shall be prepared by the recon prepare, e.g. the reconnode_backprojectionprepare. This
% function is for the cases without BP nodes for only postprocess pipelines, but these values is not so accurate as which
% prepared by the recon preapres.

if isempty(cimage)
    cimage = struct();
end

cimage.imagesize = protocol.imagesize;
cimage.center = protocol.reconcenter;
cimage.imageincrement = protocol.imageincrement;
cimage.imagethickness = protocol.imagethickness;
cimage.reconmethod = protocol.scan;
cimage.voxelsize = protocol.reconFOV/min(protocol.imagesize);
cimage.delta_radius = 1;
cimage.effFOV = protocol.effectiveFOV;
if strcmp(protocol.rawdatastyle, 'complex')
    cimage.databasis = 2;
else
    cimage.databasis = 1;
end

end