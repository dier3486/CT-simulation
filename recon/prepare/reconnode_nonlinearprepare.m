function [dataflow, prmflow, status] = reconnode_nonlinearprepare(dataflow, prmflow, status)
% corr node, non-linear correction prepare
% [dataflow, prmflow, status] = reconnode_nonlinearprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% calibration table
nonlcorr = prmflow.corrtable.(nodename);
% check focal number
if isfield(prmflow.raw, 'Nfocal')
    Nfocal = prmflow.raw.Nfocal;
else
    Nfocal = 1;
end
if Nfocal~=nonlcorr.focalnumber
    status.errormsg = sprintf('The focal spots'' number is not matching with the %s calibration table!', nodename);
    status.errorcode = 221;
    % force to use an un-matching calibration table.
end

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
    prmflow.corrtable.(nodename) = putinGPU(prmflow.corrtable.(nodename));
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end