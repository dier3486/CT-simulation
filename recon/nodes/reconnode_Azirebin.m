function [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status)
% recon node, Azi rebin for axial
% [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status);
% fly-focal is not supported yet

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nshot = prmflow.recon.Nshot;
% I know Nview=Nviewprot*Nshot, for axial

% rebin is prepared
rebin = prmflow.rebin;

% Azi rebin
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
for ishot = 1:Nshot
    A = zeros(Npixel*Nslice, Nviewprot, 'single');
    viewindex = (1:Nviewprot) + (ishot-1)*Nviewprot;
    A(rebin.vindex1_azi) = ...
        dataflow.rawdata(:, viewindex).*repmat(1-rebin.interalpha_azi, 1, Nviewprot);
    A(rebin.vindex2_azi) = A(rebin.vindex2_azi) + ...
        dataflow.rawdata(:, viewindex).*repmat(rebin.interalpha_azi, 1, Nviewprot);
    % start angle for first rebin view, replace the rawdata
    dataflow.rawdata(:, viewindex) = ...
        [A(:, (rebin.startvindex : Nviewprot)) A(:, (1 : rebin.startvindex-1))];
end
% viewangle
viewangle = rehsape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
startviewangle = viewangle(rebin.startvindex, :);
dataflow.rawhead.viewangle = [viewangle(rebin.startvindex : Nviewprot, :); viewangle(1 : rebin.startvindex-1, :)];

% prm to return
prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = rebin.delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end