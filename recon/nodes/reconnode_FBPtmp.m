function [dataflow, prmflow, status] = reconnode_FBPtmp(dataflow, prmflow, status)
% recon node, FBP 
% [dataflow, prmflow, status] = reconnode_FBPtmp(dataflow, prmflow, status);
% hard code for temporay using

% BP parameter
recon = prmflow.recon;

% hardcode
recon.N = 512;
recon.FOV = 500;
recon.h = recon.FOV/recon.N;

% read filter
% fid = fopen('D:\matlab\data\simulation\BodySoft_QDO.res');
fid = fopen('E:\data\simulation\kernel\BodySoft.bin.res');
myfilter = fread(fid, inf, 'single=>single');
fclose(fid);

dataflow.rawdata = reshape(dataflow.rawdata, recon.Npixel, recon.Nslice, recon.Nviewprot, recon.Nshot);

dataflow.image = zeros(recon.N, recon.N, recon.Nslice*recon.Nshot, 'single');
for ishot = 1:recon.Nshot
    recon.viewangle = (0:recon.Nviewprot-1).*recon.delta_view + recon.startviewangle(ishot);
    sliceindex = (1:recon.Nslice) + (ishot-1).*recon.Nslice;
    dataflow.image(:,:, sliceindex) = filterbackproj2D(squeeze(dataflow.rawdata(:, :, :, ishot)), recon, myfilter);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end