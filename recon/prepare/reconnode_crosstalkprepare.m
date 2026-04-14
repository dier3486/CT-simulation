function [dataflow, prmflow, status] = reconnode_crosstalkprepare(dataflow, prmflow, status)
% corr node, cross-tale correction prepare
% [dataflow, prmflow, status] = reconnode_crosstalkprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

if ~isfield(prmflow.pipe.(nodename), 'weight')
    prmflow.pipe.(nodename).weight = single(1.0);
end
if ~isfield(prmflow.pipe.(nodename), 'istointensity')
    % default, do crosstalk in intensity
    prmflow.pipe.(nodename).istointensity = true;
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % viewcommonfactor
    prmflow.pipe.(nodename).pipeline.viewcommonfactor = prmflow.raw.Nfocal;
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

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    % put corrtable to GPU
    [prmflow.corrtable.(nodename), prmflow.pipe.(nodename).weight] = ...
        putinGPU(prmflow.corrtable.(nodename), prmflow.pipe.(nodename).weight);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end