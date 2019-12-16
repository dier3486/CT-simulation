function [dataflow, prmflow, status] = reconnode_Axialrebin_backup(dataflow, prmflow, status)
% recon node, Axial rebin, back up
% [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nshot = prmflow.recon.Nshot;
% I know Nview=Nviewprot*Nshot, for axial
focalspot = prmflow.system.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.system.Nfocal;
% fly-focal is not supported yet
rebinpipe = prmflow.pipe.(status.nodename);

% isDQO
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
rebin = rebinprepare(detector, focalposition, Nviewprot, isQDO);

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

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nview);

% QDO reorder
if isQDO
    A_QDO = zeros(rebin.Npixel, Nslice*Nview/2, 'single');
    % I know in QDO the rebin.Npixel~=Npixel
    Savail = ~isnan(rebin.QDOorder);
    % QDO phase1
    index_s1 = (1:Nslice*Nviewprot/2)' + (0:Nshot-1).*Nslice*Nviewprot;
    A_QDO(rebin.QDOorder(Savail(:,1), 1), :) = dataflow.rawdata(Savail(:,1), index_s1(:));
    % QDO phase2
    index_s2 = (Nslice*Nviewprot/2+1:Nslice*Nviewprot)' + (0:Nshot-1).*Nslice*Nviewprot;
    A_QDO(rebin.QDOorder(Savail(:,2), 2), :) = dataflow.rawdata(Savail(:,2), index_s2(:));
    % replace the rawdata
    dataflow.rawdata = A_QDO;
end

% radial rebin
dataflow.rawdata = dataflow.rawdata(rebin.radialindex, :).*(1-rebin.interalpha_rad) + ...
    dataflow.rawdata(rebin.radialindex+1,:).*rebin.interalpha_rad;

% to return
prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = rebin.Nviewprot;
prmflow.recon.Nview = rebin.Nviewprot*Nshot;
prmflow.recon.midchannel = rebin.midchannel;
prmflow.recon.delta_d = rebin.delta_d;
prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = rebin.delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end