function [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status)
% recon node, BP (after filter)
% [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);

% BP prepare
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

% BP parameters
BPprm = prmflow.pipe.(status.nodename);
if isfield(BPprm, 'interp')
    interp = BPprm.interp;
else
    interp = 'linear';
end

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end

% recon parameters
Nshot = prmflow.recon.Nshot;
Nviewprot = prmflow.recon.Nviewprot;
delta_view = prmflow.recon.delta_view;
startviewangle = double(prmflow.recon.startviewangle);
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d;
FOV = prmflow.recon.FOV;
N = prmflow.recon.imagesize;
hond = FOV/N/delta_d;
imagecenter = prmflow.recon.imagecenter./delta_d;
maxr = prmflow.recon.maxradius;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);
% ini image
dataflow.image = zeros(N, N, Nslice*Nshot, 'single');
for ishot = 1 : Nshot 
    % view angle
    theta = (0:Nviewprot-1).*delta_view + startviewangle(ishot) + pi/2;
    % theta = viewangle+pi/2, which is the polar angle of the beams on x-y plane, compatible with matlab iradon.m
    % NOTE: theta shall in double to call iradon
    sliceindex = (1:Nslice) + (ishot-1).*Nslice;
    if GPUonoff
        dataflow.image(:,:, sliceindex) = backproj2D_GPU(dataflow.rawdata(:, :, :, ishot), theta, midchannel, hond, N, ...
            imagecenter(sliceindex, :), maxr);
    else
        dataflow.image(:,:, sliceindex) = backproj2D_2(dataflow.rawdata(:, :, :, ishot), theta, midchannel, hond, N, interp, ...
            imagecenter(sliceindex, :));
    end
end

reorderflag = prmflow.protocol.couchdirection < 0;
dataflow.image = imagereorder(dataflow.image, Nslice, reorderflag);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end