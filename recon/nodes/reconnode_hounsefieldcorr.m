function [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status)
% recon node, Hounsefield Units correction
% [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status);

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
% parameters set in pipe
HCprm = prmflow.pipe.(status.nodename);

if isfield(HCprm, 'HCscale')
    HCscale = HCprm.HCscale;
else
    HCscale = 1000;
end

% scale
dataflow.rawdata = dataflow.rawdata.*HCscale;

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    HUcorr = prmflow.corrtable.(status.nodename);
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview).*HUcorr.main(:)';
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end