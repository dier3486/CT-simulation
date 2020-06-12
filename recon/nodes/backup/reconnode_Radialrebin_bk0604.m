function [dataflow, prmflow, status] = reconnode_Radialrebin(dataflow, prmflow, status)
% recon node, Axial rebin 
% [dataflow, prmflow, status] = reconnode_Radialrebin(dataflow, prmflow, status);

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nshot = prmflow.recon.Nshot;
% I know Nview=Nviewprot*Nshot, for axial

% rebin is prepared
rebin = prmflow.rebin;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nview);

% QDO reorder
if rebin.isQDO
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

% prm to return
prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = rebin.Nviewprot;
prmflow.recon.Nview = rebin.Nviewprot*Nshot;
prmflow.recon.midchannel = rebin.midchannel;
prmflow.recon.delta_d = rebin.delta_d;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end