function [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status)
% recon node, beamharden correction
% [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status);

% parameters to use in prmflow
Nview = prmflow.recon.Nview;
% Npixel = prmflow.recon.Npixel;
% Nslice = prmflow.recon.Nslice;

% calibration table
bhcorr = prmflow.corrtable.(status.nodename);
bhorder = bhcorr.order;
bhpoly = reshape(bhcorr.main, [], bhorder);

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
dataflow.rawdata = iterpolyval(bhpoly, dataflow.rawdata);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end