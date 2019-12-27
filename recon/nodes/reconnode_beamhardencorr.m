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

% % test nonlinear
% p30 = load('E:\data\rawdata\bhtest\p30.mat');
% p = p30.p;
% p(p(:,2)==0, 2) = 1.0;
% dataflow.rawdata = reshape(dataflow.rawdata, size(p,1), []);
% dataflow.rawdata = dataflow.rawdata.^2.*p(:,1) + dataflow.rawdata.*p(:,2);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end