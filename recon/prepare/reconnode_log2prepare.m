function [dataflow, prmflow, status] = reconnode_log2prepare(dataflow, prmflow, status)
% prepare node, prepare of log2
% [dataflow, prmflow, status] = reconnode_log2prepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

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
end

end