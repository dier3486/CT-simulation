function [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status)
% corr node, air correction prepare
% [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% calibration table
aircorr = prmflow.corrtable.(nodename);
% check focal number
if isfield(prmflow.raw, 'Nfocal')
    Nfocal = prmflow.raw.Nfocal;
else
    Nfocal = 1;
end
if Nfocal~=aircorr.focalnumber
    status.errormsg = 'The focal spots'' number is not matching with the air calibration table!';
    status.errorcode = 221;
    % force to use an un-matching calibration table.
end

Nref = aircorr.refnumber;
% reference error cut (ref-error-cut)
if isfield(aircorr, 'referrcut')
    aircorr.referrcut = reshape(aircorr.referrcut, Nref, aircorr.focalnumber);
    if Nfocal==aircorr.focalnumber
        referrcut = aircorr.referrcut;
    else
        referrcut = zeros(Nref, Nfocal, 'single');
        for ifocal = 1:Nfocal
            referrcut(:, ifocal) = aircorr.referrcut(:, mod(ifocal-1, aircorr.focalnumber) + 1);
        end
    end
else
    referrcut = zeros(Nref, Nfocal, 'single');
end

% save the prepared parameters
prmflow.correction.air = struct();
prmflow.correction.air.referrcut = single(referrcut);
% refpixel, refnumber
prmflow.correction.air.refpixel = single(aircorr.refpixel);
prmflow.correction.air.refnumber = single(aircorr.refnumber);
% I know the referrcut and other ref* will be useding reference correction but (could) not in air correction.

% airKVmA
prmflow.correction.air.airKVmA = aircorr.referenceKVmA;

% pipe line
if pipeline_onoff
    % pipeline console paramters, default
    prmflow.pipe.(nodename).pipeline = struct();
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end