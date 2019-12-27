function [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status)
% recon node, filter design and conv
% [dataflow, prmflow, status] = reconnode_filter(dataflow, prmflow, status);

% prm
Npixel = prmflow.recon.Npixel;
delta_d = prmflow.recon.delta_d;

% filter
filter = prmflow.pipe.(status.nodename);

% design filter
prmflow.recon.filter = loadfilter(filter, Npixel, delta_d);

% conv
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, []);
% fill zero
dataflow.rawdata(length(prmflow.recon.filter), :) = 0;
% fft
dataflow.rawdata = fft(dataflow.rawdata);
% time
dataflow.rawdata = dataflow.rawdata.*prmflow.recon.filter;
% ifft
dataflow.rawdata = ifft(dataflow.rawdata, 'symmetric');
% kick filled zero
dataflow.rawdata(Npixel+1:end,:) = [];
% done

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end