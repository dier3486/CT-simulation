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
nonlpoly = reshape(nonlcorr.main, [], nonlorder);

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
dataflow.rawdata = iterpolyval(nonlpoly, dataflow.rawdata);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end