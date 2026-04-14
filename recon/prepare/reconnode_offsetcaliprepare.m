function [dataflow, prmflow, status] = reconnode_offsetcaliprepare(dataflow, prmflow, status)
% prepare node, prepare of offset cali
% [dataflow, prmflow, status] = reconnode_offsetcaliprepare(dataflow, prmflow, status);

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
end

% initial the buffer
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
dataflow.buffer.(nodename).Nview = 0;
dataflow.buffer.(nodename).miu = zeros(Npixel*Nslice,1, 'single');
dataflow.buffer.(nodename).variance = zeros(Npixel*Nslice,1, 'single');

end