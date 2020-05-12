function [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status)
% recon node, nonlinear correction
% [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status);
% the algorithm is exactly same as beamharden corr, so we suggest to assign beamharden and nonlinear correction in one node in
% nodesentry.m

% parameters to use in prmflow
Nview = prmflow.recon.Nview;

% calibration table
nonlcorr = prmflow.corrtable.(status.nodename);
nonlorder = nonlcorr.order;
% DFS
if isfield(nonlcorr, 'focalnumber') && nonlcorr.focalnumber
    Nfocal = nonlcorr.focalnumber;
else
    Nfocal = 1;
end
nonlpoly = reshape(nonlcorr.main, [], nonlorder, Nfocal);

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
for ifocal = 1:Nfocal
    dataflow.rawdata(:, ifocal:Nfocal:end) = iterpolyval(nonlpoly(:, :, ifocal), dataflow.rawdata(:, ifocal:Nfocal:end));
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end