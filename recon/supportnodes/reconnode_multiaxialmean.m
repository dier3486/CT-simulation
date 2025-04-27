function [dataflow, prmflow, status] = reconnode_multiaxialmean(dataflow, prmflow, status)
% mean multi-rotation of rawdata
% [dataflow, prmflow, status] = reconnode_multiaxialmean(dataflow, prmflow, status);

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_multiaxialmeanprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

Nshot = prmflow.raw.Nshot;
Nviewprot = prmflow.raw.Nviewprot;
Nmulti = prmflow.raw.Nmulti;

% mean
dataflow.rawdata = reshape((mean(reshape(dataflow.rawdata, [], Nviewprot, Nmulti, Nshot), 3)), [], Nviewprot*Nshot);

headfields = fieldnames(dataflow.rawhead);
for ii = 1 : length(headfields)
    switch headfields{ii}
        case 'Reading_Number'
            tmp = reshape(dataflow.rawhead.(headfields{ii}), Nviewprot, Nmulti, Nshot);
            dataflow.rawhead.(headfields{ii}) = reshape(tmp(:,1,:), 1, Nviewprot*Nshot);
        otherwise
            dataflow.rawhead.(headfields{ii}) = reshape( ...
                mean(reshape(dataflow.rawhead.(headfields{ii}), [], Nviewprot, Nmulti, Nshot), 3), [], Nviewprot*Nshot);
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end