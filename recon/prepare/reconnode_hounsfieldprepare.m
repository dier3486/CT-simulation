function [dataflow, prmflow, status] = reconnode_hounsfieldprepare(dataflow, prmflow, status)
% corr node, Hounsfield correction prepare
% [dataflow, prmflow, status] = reconnode_hounsfieldprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end