function [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status)
% recon node, filter design and conv
% [dataflow, prmflow, status] = reconnode_filter(dataflow, prmflow, status);

% prm
Npixel = prmflow.recon.Npixel;
delta_d = prmflow.recon.delta_d;

% filter
filter = prmflow.pipe.(status.nodename);

% design filter
if isfield(filter, 'name')  || ifchar(filter)
    % design filter by filter's name
    if isfield(filter, 'freqscale')
        freqscale = filter.freqscale;
    else
        freqscale = 1.0;
    end
    if ischar(filter)
        filtname = filter;
    else
        filtname = filter.name;
    end
    H = filterdesign(filtname, Npixel, delta_d, freqscale);
elseif isfield(filter, 'file')
    % load filter from a file
    if ~exist(filter.file, 'file')
        error('Can not find filter file: %s', filter.file);
    end
    fid = fopen(filter.file, 'r');
    H = fread(fid, inf, 'single');
    fclose(fid);
else
    error('Filter is not defined!');
end
prmflow.recon.filter = H;

% conv
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, []);
% fill zero
dataflow.rawdata(length(H), :) = 0;
% fft
dataflow.rawdata = fft(dataflow.rawdata);
% time
dataflow.rawdata = dataflow.rawdata.*H;
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