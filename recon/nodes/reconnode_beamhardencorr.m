function [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status)
% recon node, beamharden correction
% [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status);

% parameters to use in prmflow
Nview = prmflow.recon.Nview;

% calibration table
bhcorr = prmflow.corrtable.(status.nodename);
bhorder = bhcorr.order;
bhpoly = reshape(bhcorr.main, [], bhorder);

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
Dchk_corr = ones(size(dataflow.rawdata), 'single');
for ii = 1:bhorder
    Dchk_corr = (Dchk_corr.*dataflow.rawdata).*bhpoly(:,ii) + 1.0;
end
dataflow.rawdata = Dchk_corr - 1.0;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end