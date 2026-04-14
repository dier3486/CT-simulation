function [dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status)
% prepare node, data buffer prepare
% [dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
if ~status.pipeline.(nodename).pipeline_onoff
    status.jobdone = false;
    status.errorcode = -1;
    status.errormsg = sprintf('The %s can be called only while pineline on.', nodename);
    return;
end

% default paramters
% if to copy the buffer back to dataflow (after previous works done)
if ~isfield(nodeprm, 'copytodataflow')
    prmflow.pipe.(nodename).copytodataflow = false;
end
% if to stuck the pipeline until previous works done, used to link a pipeline node with a non-pipeline node.
if ~isfield(nodeprm, 'stuck')
    prmflow.pipe.(nodename).stuck = false;
end
% never carry the databuffer
prmflow.pipe.(nodename).pipeline.iscarried = false;

% go to default prepare
[dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status);
% skip the default prepare in reconnode_pipelineprepare
status.currentjob.prepared = true;

% the data to be buffered
if isfield(nodeprm, 'alldata') && nodeprm.alldata
    bufferfields = dataflow.pipepool.(nodename)(1).datafields;
else
    if isfield(nodeprm, 'bufferfields')
        bufferfields = regexp(char(nodeprm.bufferfields), '(, +)|(,)', 'split');
    else
        bufferfields = dataflow.pipepool.(nodename)(1).datafields;
    end
end
% exclude fields
if isfield(nodeprm, 'excludefields')
    excludefields = regexp(char(nodeprm.excludefields), '(, +)|(,)', 'split');
else
    excludefields = {};
end
bufferfields = setdiff(bufferfields, excludefields);

% ini buffer
Npool = length(dataflow.pipepool.(nodename));
dataflow.buffer.(nodename).pool(1 : Npool) = status.defaultpool;
% buffer's datafields is bufferfields
dataflow.buffer.(nodename).pool(1).datafields = bufferfields;
dataflow.buffer.(nodename).pool(1).databasis = dataflow.pipepool.(nodename)(1).databasis;
% copy pipepool.data to buffer.data
dataflow.buffer.(nodename).pool(1).data = struct();
for ifield = bufferfields
    if isfield(dataflow.pipepool.(nodename).data, ifield{1})
        if ~isa(dataflow.pipepool.(nodename)(1).data.(ifield{1}), 'lib.pointer')
            dataflow.buffer.(nodename).pool(1).data.(ifield{1}) = gather(dataflow.pipepool.(nodename)(1).data.(ifield{1}));
        else
            % to gather the lib.pointers
            datasize = dataflow.pipepool.(nodename)(1).data.([ifield{1}, '_datasize']);
            databasis = dataflow.pipepool.(nodename)(1).data.([ifield{1}, '_databasis']);
            if databasis ~= 2
                dataflow.buffer.(nodename).pool(1).data.(ifield{1}) = zeros(datasize*databasis, 0, 'single');
            else
                dataflow.buffer.(nodename).pool(1).data.(ifield{1}) = complex(zeros(datasize, 0, 'single'));
            end
            % only single now
        end
    else
        dataflow.buffer.(nodename).pool(1).data.(ifield{1}) = [];
    end
end
% if Npool > 1
for ii = 2:Npool
    % mirror the datafields and databasis
    dataflow.buffer.(nodename).pool(ii).datafields = dataflow.pipepool.(nodename)(ii).datafields;
    dataflow.buffer.(nodename).pool(ii).databasis = dataflow.pipepool.(nodename)(ii).databasis;
    % copy (gather) the other pools
    for ifield = dataflow.pipepool.(nodename)(ii).datafields
        dataflow.buffer.(nodename).pool(ii).data.(ifield{1}) = gather(dataflow.pipepool.(nodename)(ii).data.(ifield{1}));
    end
end

% GPU on/off
% default GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
    prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
end

% carried? can not be carried yet
% if prmflow.protocol.tocarrythepools     % default is true
%     prmflow.pipe.(nodename).pipeline.iscarried = true;
%     % default was false
% end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end