function [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status)
% recon node, crosstalk correction
% [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status);

% parameters to use in prmflow
Nview = prmflow.recon.Nview;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;

% parameters set in pipe
crossprm = prmflow.pipe.(status.nodename);
if isfield(crossprm, 'weight')
    weight = crossprm.weight;
else
    weight = 1.0;
end

% calibration table
crscorr = prmflow.corrtable.(status.nodename);
crsorder = crscorr.order;
crsval = reshape(crscorr.main, [], bhorder);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
% to intensity
dataflow.rawdata = 2.^(-dataflow.rawdata);
% correct
if crsorder == 1
    crsval = reshape(crsval, Npixel, Nslice);
    dataflow.rawdata(2:end-1,:,:) = dataflow.rawdata(2:end-1,:,:)+(dataflow.rawdata(3:end,:,:)-dataflow.rawdata(1:end-2,:,:)).*crsval(2:end-1,:).*weight;
else
    error('Currently thhe crosstalk correction only support 1-order method.');
end
% min cut
minval = 2^-32;
dataflow.rawdata(dataflow.rawdata<minval) = minval;
% log2
dataflow.rawdata = -log2(dataflow.rawdata);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end