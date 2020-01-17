function [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status)
% recon node, filter design and conv
% [dataflow, prmflow, status] = reconnode_filter(dataflow, prmflow, status);

% prm
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
delta_d = prmflow.recon.delta_d;

% filter
filter = prmflow.pipe.(status.nodename);

% design filter
prmflow.recon.filter = loadfilter(filter, Npixel, delta_d);
Hlen = length(prmflow.recon.filter);

% fill
if isfield(filter, 'fillup') && filter.fillup
    % fill up
    if isfield(dataflow.rawhead, 'refblock')
        blkvindex = any(dataflow.rawhead.refblock, 1);
    else
        blkvindex = [];
    end
    mid_u = prmflow.recon.midchannel;
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
    A = zeros(Hlen, Nslice, Nview);
    for ii = 1:Nslice
        [A(:, ii, :), n_left] = translatefillup(squeeze(dataflow.rawdata(:, ii, :)), Hlen, mid_u, blkvindex);
    end
    dataflow.rawdata = reshape(A, Hlen, []);
else
    % fill zero
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, []);
    dataflow.rawdata(Hlen, :) = 0;
    n_left = 0;
end

% conv
% fft
dataflow.rawdata = fft(dataflow.rawdata);
% time
dataflow.rawdata = dataflow.rawdata.*prmflow.recon.filter;
% ifft
dataflow.rawdata = ifft(dataflow.rawdata, 'symmetric');
% kick filled zero
% dataflow.rawdata(Npixel+1:end,:) = [];
dataflow.rawdata = dataflow.rawdata((1:Npixel)+n_left, :);
% done

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end