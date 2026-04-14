function [dataflow, prmflow, status] = reconnode_smartairprepare(dataflow, prmflow, status)
% corr node, air correction prepare
% [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% calibration table
smartaircorr = prmflow.corrtable.(nodename);
% check focal number
if isfield(prmflow.raw, 'Nfocal')
    Nfocal = prmflow.raw.Nfocal;
else
    Nfocal = 1;
end
if Nfocal~=smartaircorr.focalnumber
    status.errormsg = 'The focal spots'' number is not matching with the air calibration table!';
    status.errorcode = 221;
    % force to use an un-matching calibration table.
end

%{
Nref = smartaircorr.refnumber;
% reference error cut (ref-error-cut)
if isfield(smartaircorr, 'referrcut')
    smartaircorr.referrcut = reshape(smartaircorr.referrcut, Nref, smartaircorr.focalnumber);
    if Nfocal==smartaircorr.focalnumber
        referrcut = smartaircorr.referrcut;
    else
        referrcut = zeros(Nref, Nfocal, 'single');
        for ifocal = 1:Nfocal
            referrcut(:, ifocal) = smartaircorr.referrcut(:, mod(ifocal-1, smartaircorr.focalnumber) + 1);
        end
    end
else
    referrcut = zeros(Nref, Nfocal, 'single');
end

% save the prepared parameters
prmflow.correction.smartair = struct();
prmflow.correction.smartair.referrcut = single(referrcut);
% refpixel, refnumber
prmflow.correction.smartair.refpixel = single(smartaircorr.refpixel);
prmflow.correction.smartair.refnumber = single(smartaircorr.refnumber);
% I know the referrcut and other ref* will be useding reference correction but (could) not in air correction.

% airKVmA
prmflow.correction.smartair.airKVmA = smartaircorr.referenceKVmA;
%}

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