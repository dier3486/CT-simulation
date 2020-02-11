function [dataflow, prmflow, status] = reconnode_datamean(dataflow, prmflow, status)
% support node, mean of rawdata
% [dataflow, prmflow, status] = reconnode_datamean(dataflow, prmflow, status);

% parameters set in pipe
meanprm = prmflow.pipe.(status.nodename);

% skipviews
if isfield(meanprm, 'viewskip') && ~isempty(meanprm.viewskip)
    viewskip = meanprm.viewskip;
else
    viewskip = 0;
end

% mean
rawsize = size(dataflow.rawdata);
dataflow.rawdata = reshape(dataflow.rawdata, [], rawsize(end));
dataflow.rawdata = mean(dataflow.rawdata(:, viewskip+1:end), 2);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end