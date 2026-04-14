function [dataflow, prmflow, status] = reconnode_aircaliprepare(dataflow, prmflow, status)
% prepare node, prepare of air cali
% [dataflow, prmflow, status] = reconnode_aircaliprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

if ~isfield(nodeprm, 'corrversion')
    prmflow.pipe.(nodename).corrversion = 'v1.0';
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % default GPU off
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 0;
    end
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end
    % prmflow.pipe.(nodename).pipeline.iscarried = false;
    prmflow.pipe.(nodename).pipeline.inputminlimit = prmflow.raw.Nviewprot;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end