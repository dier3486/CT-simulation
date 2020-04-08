function [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status)
% recon node, BP (after filter)
% [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);

% BP parameters
BPprm = prmflow.pipe.(status.nodename);
if isfield(BPprm, 'FOV')
    FOV = BPprm.FOV;
else
    FOV = 500;
end
if isfield(BPprm, 'imagesize')
    N = BPprm.imagesize;
else
    N = 512;
end
if isfield(BPprm, 'interp')
    interp = BPprm.interp;
else
    interp = 'linear';
end

% recon parameters
Nshot = prmflow.recon.Nshot;
Nviewprot = prmflow.recon.Nviewprot;
delta_view = prmflow.recon.delta_view;
startviewangle = prmflow.recon.startviewangle;
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d;
hond = FOV/N/delta_d;

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
    dataflow.image(:,:, sliceindex) = backproj2D_1(dataflow.rawdata(:, :, :, ishot), theta, midchannel, hond, N, interp);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end