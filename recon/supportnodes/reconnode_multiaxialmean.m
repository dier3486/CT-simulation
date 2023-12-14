function [dataflow, prmflow, status] = reconnode_multiaxialmean(dataflow, prmflow, status)
% mean multi-rotation of rawdata
% [dataflow, prmflow, status] = reconnode_multiaxialmean(dataflow, prmflow, status);

% parameters set in pipe
nodeprm = prmflow.pipe.(status.nodename);

% if to mean all the shots
if isfield(nodeprm, 'meanshots')
    meanshots = nodeprm.meanshots;
else
    meanshots = false;
end

% parameters to use in prmflow
if meanshots
    Nshot = 1;
else
    Nshot = prmflow.raw.Nshot;
end
Nview = prmflow.raw.Nview;
Nviewprot = prmflow.raw.Nviewprot;
Nmulti = Nview/Nviewprot;

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

% prmflow.raw after mean ??
prmflow.raw.Nview = Nviewprot * Nshot;
prmflow.raw.Nshot = Nshot;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end